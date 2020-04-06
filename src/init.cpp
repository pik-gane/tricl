/*
 * init.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include <assert.h>
#include <iostream>

#include "debugging.h"
#include "data_model.h"
#include "global_variables.h"
#include "probability.h"
#include "entity.h"
#include "link.h"
#include "event.h"
#include "gexf.h"
#include "init.h"
#include "io.h"

using namespace std;

int n_rats = 0;

// arrays of parameters:
rate _inflt2attempt_rate[MAX_N_INFLT];
probunit _inflt2delta_probunit[MAX_N_INFLT];
entity_type _e2et[MAX_N_E];

// derived constants:
entity max_e = -1; // largest entity id in use
set<entity> es;
unordered_map<entity_type_pair, set<relationship_or_action_type>> ets2relations;  // possible relations

// variable data:

timepoint current_t = 0;
event current_ev = {};
event_data* current_evd_ = NULL;

// network state:
unordered_map<entity_type, vector<entity>> et2es = {};  // kept to equal inverse of v2vt
unordered_map<entity, outleg_set> e2outs = {};
unordered_map<entity, inleg_set> e2ins = {};
int n_links = 0; // total no. of current (non-id.) links incl. inverse relationships

// event data:
unordered_map<event, event_data> ev2data = {};
map<timepoint, event> t2be = {};  // kept to equal inverse of ev2data.t


void init_data ()
{
    if (debug) {
        verbose = true;
    }
    if (verbose) {
        quiet = false;
    }
    // init with zeroes:
    for (int i=0; i<MAX_N_INFLT; i++) {
        _inflt2delta_probunit[i] = _inflt2attempt_rate[i] = 0;
    }
    // store actual values:
    for (auto& [inflt, ar] : inflt2attempt_rate) {
        assert (!( inflt.evt.ec == EC_EST && ( inflt.at.rat12 == NO_RAT || inflt.at.rat23 == NO_RAT ) && !(inflt.at == NO_ANGLE)));
        _inflt2attempt_rate[INFLT(inflt)] = ar;
    }
    for (auto& [inflt, spu] : inflt2delta_probunit) {
        assert (!( inflt.evt.ec == EC_EST && ( inflt.at.rat12 == NO_RAT || inflt.at.rat23 == NO_RAT ) ));
        _inflt2delta_probunit[INFLT(inflt)] = spu;
    }
    for (auto& [e, et] : e2et) {
        _e2et[E(e)] = et;
    }
}

void init_entities ()
{
    // inspect pre-registered entities:
    auto vt2remaining_n = et2n;
    for (auto& [e, et] : e2et) {
        assert(e >= 0);
        es.insert(e);
        if (e > max_e) max_e = e;
        if (e2label.count(e) == 0) e2label[e] = to_string(e);
        if (e2outs.count(e) == 0) e2outs[e] = {};
        if (e2ins.count(e) == 0) e2ins[e] = {};
        if (et2es.count(et) == 0) et2es[et] = {};
        _e2et[E(e)] = et;
        et2es[et].push_back(e);
        vt2remaining_n[et]--;
    }

    // generate remaining entities:
    for (auto& [et, n] : vt2remaining_n) {
        cout << " entity type \"" << et2label[et] << "\" has " << et2n[et] << " entities" << endl;
        if (n < et2n[et]) {
            cout << "  of which were preregistered:" << endl;
            for (auto& e : et2es[et]) cout << "   " << e2label[e] << endl;
        }
        assert (n >= 0);
        while (n > 0) {
            add_entity(et);
            n--;
        }
        assert ((entity) et2es[et].size() == et2n[et]);
    }
}

void init_relationship_or_action_types ()
{
    // verify symmetry of relationship inversion map:
    cout << " relationship types (with inverses) and action types:" << endl;
    assert (rat2inv[RT_ID] == RT_ID);
    for (auto& [r, la] : rat2label) {
        auto inv = rat2inv[r];
        if (inv == NO_RAT) {
            if (r_is_action_type[r]) {
                cout << "  non-symmetric action type \"" << la << "\" (no inverse)" << endl;
            } else {
                cout << "  non-symmetric relationship type \"" << la << "\" (no inverse)" << endl;
            }
        } else {
            if (inv == r) {
                if (r_is_action_type[r]) {
                    cout << "  symmetric action type \"" << la << "\"" << endl;
                } else {
                    cout << "  symmetric relationship type \"" << la << "\"" << endl;
                }
            } else {
                assert (rat2inv[inv] == r);
                assert (r_is_action_type[r] == r_is_action_type[inv]);
                if (r_is_action_type[r]) {
                    cout << "  non-symmetric action type \"" << la << "\" (inverse: \"" << rat2label[inv] << "\")" << endl;
                } else {
                    cout << "  non-symmetric relationship type \"" << la << "\" (inverse: \"" << rat2label[inv] << "\")" << endl;
                }
            }
        }
    }

    // register possible relationship types by entity type pair, and compute probunits:
    cout << " base success probabilities:" << endl;
    for (auto& [evt, pu] : evt2base_probunit) {
        auto rat13 = evt.rat13;
        assert (rat13 != RT_ID);
        entity_type_pair ets = { evt.et1, evt.et3 };
        if (ets2relations.count(ets) == 0) ets2relations[ets] = {};
        ets2relations[ets].insert(rat13);
        cout << "  " << evt << ": " << probunit2probability(pu, evt2left_tail.at(evt), evt2right_tail.at(evt)) << endl;
    }

    n_rats = rat2label.size();
}

void init_summary_events ()
{
    cout << " initial scheduling of summary events" << endl;
    // summary events for purely spontaneous establishment without angles:
    for (auto& [ets, relations] : ets2relations) {
        auto et1 = ets.et1, et3 = ets.et3;
        for (auto& rat13 : relations) {
            event ev = { .ec = EC_EST,
                    .e1 = (entity)-et1, // in spontaneous events, fields e1 and e3 are used to store entity types with negative sign
                    .rat13 = rat13,
                    .e3 = (entity)-et3 };
            event_type evt = { .ec = EC_EST, .et1 = et1, .rat13 = rat13, .et3 = et3 };
            influence_type inflt = { .evt = evt, .at = NO_ANGLE };
            auto ar = _inflt2attempt_rate[INFLT(inflt)];
            if (ar > 0) {
                if (verbose) cout << "  " << et2label[et1] << " " << rat2label[rat13] << " " << et2label[et3] << endl;
                ev2data[ev] = { .n_angles = 0,
                        .attempt_rate = ar * et2n[et1] * et2n[et3],
                        .success_probunits = evt2base_probunit[evt] + _inflt2delta_probunit[INFLT(inflt)] };
                schedule_event(ev, &ev2data[ev], evt2left_tail.at(evt), evt2right_tail.at(evt));
            }
        }
    }
}

void init_links ()
{
    cout << " performing events that add initial links" << endl;

    // identity relationship:
    for (auto& e : es) {
        e2outs[e].insert({ RT_ID, e });
        e2ins[e].insert({ e, RT_ID });
    }

    // preregistered links:
    for (auto l : initial_links) {
        assert (!link_exists(l));
        assert (l.rat13 != RT_ID);
        event ev = { .ec = EC_EST, .e1 = l.e1, .rat13 = l.rat13, .e3 = l.e3 };
        conditionally_remove_event(ev);
        perform_event(ev);
    }

    // random links:

    // block model:
    unordered_map<entity, int> e2block = { };
    for (auto& e : es) {
        e2block[e] = floor(uniform(random_variable) * et2n_blocks[_e2et[E(e)]]);
    }
    for (auto& [lt, pw] : lt2initial_prob_within) { // TODO: add lt2initial_prob_between
        if (pw > 0) {
            assert (lt.rat13 != RT_ID);
            auto et1 = lt.et1, et3 = lt.et3; auto rat13 = lt.rat13;
            for (auto& e1 : et2es[et1]) {
                for (auto& e3 : et2es[et3]) {
                    if ((e3 != e1) && ((e3 > e1) || (rat2inv[rat13] != rat13))) {
                        float p = (e2block[e1] == e2block[e3]) ? lt2initial_prob_within[lt] : lt2initial_prob_between[lt];
                        do_random_link(p, e1, rat13, e3);
                    }
                }
            }
        }
    }

    // random geometric model:
    unordered_map<entity, vector<float>> e2coords = { };
    for (auto& [et, dim] : et2dim) {
        for (auto& e : et2es[et]) {
            for (int d=0; d<dim; d++) {
                e2coords[e].push_back(uniform(random_variable));
            }
        }
    }
    for (auto& [lt, ex] : lt2spatial_decay) {
        assert (lt.rat13 != RT_ID);
        auto et1 = lt.et1, et3 = lt.et3; auto rat13 = lt.rat13; auto dim = et2dim[et1];
        assert (et2dim[et3] == dim);
        for (auto& e1 : et2es[et1]) {
            for (auto& e3 : et2es[et3]) {
                if ((e3 != e1) && ((e3 > e1) || (rat2inv[rat13] != rat13))) {
                    float dist2 = 0.0;
                    for (int d=0; d<dim; d++) {
                        dist2 += pow(e2coords[e1][d] - e2coords[e3][d], 2);
                    }
                    float p = exp(-ex * sqrt(dist2));
                    do_random_link(p, e1, rat13, e3);
                }
            }
        }
    }
    if (debug) verify_angle_consistency();
}

void init ()
{
    cout << "INITIALIZING..." << endl;
    init_data();
    init_entities();
    init_relationship_or_action_types();
    init_summary_events();
    init_links();
    init_gexf();
    if (debug) {
        dump_data();
        verify_data_consistency();
    }
    cout << "...INITIALIZATION FINISHED." << endl << endl;
}


