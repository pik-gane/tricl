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
#include "graphviz.h"
#include "gexf.h"
#include "init.h"
#include "io.h"

using namespace std;

// parameters:
int n_rats = 0; // total no. of rats
rate _inflt2attempt_rate[MAX_N_INFLT];
probunits _inflt2delta_probunit[MAX_N_INFLT];
entity_type _e2et[MAX_N_E];

// derived constants:
set<entity> es;
unordered_map<entity_type_pair, set<relationship_or_action_type>> ets2relations;  // possible relations
unordered_map<event, rate> ev2max_sp; // max possible effective rate of summary events

// variable data:

timepoint current_t = 0;
long int n_events = 0;
event current_ev = {};
event_data* current_evd_ = NULL;

// network state:
unordered_map<entity_type, vector<entity>> et2es = {};  // kept to equal inverse of e2et
unordered_map<entity, outleg_set> e2outs = {};
unordered_map<entity, inleg_set> e2ins = {};
unordered_map<link_type, long int> lt2n = {};
long int n_links = 0, n_angles = 0;

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
    for (auto& [inflt, spu] : inflt2delta_probunits) {
        assert (!( inflt.evt.ec == EC_EST && ( inflt.at.rat12 == NO_RAT || inflt.at.rat23 == NO_RAT ) ));
        _inflt2delta_probunit[INFLT(inflt)] = spu;
    }
}

void init_entities ()
{
    // inspect pre-registered entities:
    auto et2remaining_n = et2n;
    for (auto& [e, et] : e2et) {
        if (et >= 1<<ET_BITS) throw "too many entity types (recompile with larger ET_BITS?)";
        assert(e >= 0);
//        es.insert(e);
//        if (e > max_e) max_e = e;
//        if (e2label.count(e) == 0) e2label[e] = to_string(e);
//        if (e2outs.count(e) == 0) e2outs[e] = {};
//        if (e2ins.count(e) == 0) e2ins[e] = {};
//        if (et2es.count(et) == 0) et2es[et] = {};
//        _e2et[E(e)] = et;
//        et2es[et].push_back(e);
        if (et2remaining_n[et] > 0) {
            et2remaining_n[et]--;
        } else {
            et2n[et]++;
        }
    }

    // generate remaining entities:
    for (auto& [et, n] : et2remaining_n) {
        cout << " entity type \"" << et2label[et] << "\" has " << et2n[et] << " entities" << endl;
//        if (n < et2n[et]) {
//            cout << "  of which were preregistered:" << endl;
//            for (auto& e : et2es[et]) cout << "   " << e2label[e] << endl;
//        }
        assert (n >= 0);
        while (n > 0) {
            add_entity(et, "");
            n--;
        }
        assert ((entity) et2es[et].size() == et2n[et]);
    }
    if (max_e >= 1<<E_BITS) throw "too many entities (recompile with larger E_BITS?)";
}

void init_relationship_or_action_types ()
{
    // verify symmetry of relationship inversion map:
    cout << " relationship types (with inverses) and action types:" << endl;
    assert (rat2inv.at(RT_ID) == RT_ID);
    for (auto& [r, la] : rat2label) {
        if (r >= 1<<RAT_BITS) throw "too many relationship or action types (recompile with larger RAT_BITS?)";
        if (rat2inv.count(r) == 0) rat2inv[r] = NO_RAT;
        auto inv = rat2inv[r];
        if (inv == NO_RAT) {
//            if (r_is_action_type[r]) {
//                cout << "  non-symmetric action type \"" << la << "\" (no inverse)" << endl;
//            } else {
//                cout << "  non-symmetric relationship type \"" << la << "\" (no inverse)" << endl;
//            }
        } else {
            if (inv == r) {
//                if (r_is_action_type[r]) {
//                    cout << "  symmetric action type \"" << la << "\"" << endl;
//                } else {
//                    cout << "  symmetric relationship type \"" << la << "\"" << endl;
//                }
            } else {
                assert (rat2inv.at(inv) == r);
                assert (r_is_action_type.at(r) == r_is_action_type.at(inv));
//                if (r_is_action_type[r]) {
//                    cout << "  non-symmetric action type \"" << la << "\" (inverse: \"" << rat2label[inv] << "\")" << endl;
//                } else {
//                    cout << "  non-symmetric relationship type \"" << la << "\" (inverse: \"" << rat2label[inv] << "\")" << endl;
//                }
            }
        }
    }

    // register possible relationship types by entity type pair, and compute probunits:
    cout << " base success probabilities:" << endl;
    for (auto& [evt, pu] : evt2base_probunits) {
        auto rat13 = evt.rat13;
        assert (rat13 != RT_ID);
        entity_type_pair ets = { evt.et1, evt.et3 };
        if (ets2relations.count(ets) == 0) ets2relations[ets] = {};
        ets2relations[ets].insert(rat13);
        lt2n[{evt.et1, rat13, evt.et3}] = 0;
        cout << "  " << evt << ": " << probunits2probability(pu, evt2left_tail.at(evt), evt2right_tail.at(evt)) << endl;
    }

    n_rats = rat2label.size();
}

void init_summary_events ()
{
    if (!quiet) cout << " perform initial scheduling of summary events..." << endl;
    // summary events for purely spontaneous establishment without angles:
    for (auto& [ets, relations] : ets2relations) {
        auto et1 = ets.et1, et3 = ets.et3;
        for (auto& rat13 : relations) {
            event summary_ev = { .ec = EC_EST,
                    .e1 = (entity)-et1, // in spontaneous events, fields e1 and e3 are used to store entity types with negative sign
                    .rat13 = rat13,
                    .e3 = (entity)-et3 };
            event_type evt = { .ec = EC_EST, .et1 = et1, .rat13 = rat13, .et3 = et3 };
            influence_type inflt = { .evt = evt, .at = NO_ANGLE };
            auto ar = _inflt2attempt_rate[INFLT(inflt)];
            if (ar > 0) {
                // compile maximal success units. if no influences can increase the success units,
                // this equals the base_probunits, otherwise it is infinite:
                probunits max_spu = evt2base_probunits.at(evt);
                for (auto& [inflt, pu] : inflt2delta_probunits) {
                    if ((inflt.evt == evt) && (pu > 0.0)) {
                        max_spu = INFINITY;
                        break;
                    }
                }
                ev2max_sp[summary_ev] = probunits2probability(max_spu, evt2left_tail.at(evt), evt2right_tail.at(evt));
                if (verbose) cout << "  " << et2label[et1] << " " << rat2label[rat13] << " " << et2label[et3] << endl;
                ev2data[summary_ev] = { .n_angles = 0,
                        .attempt_rate = ar * et2n[et1] * et2n[et3],
                        .success_probunits = evt2base_probunits[evt] + _inflt2delta_probunit[INFLT(inflt)] };
                schedule_event(summary_ev, &ev2data[summary_ev], evt2left_tail.at(evt), evt2right_tail.at(evt));
            }
        }
    }
    if (!quiet) cout << "  ...done."<< endl;
}

void init_links ()
{
    if (!quiet) cout << " perform events that add initial links..." << endl;

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
        perform_event(ev);  // this automatically also adds the inverse link, if any.
    }

    // random links:

    // block model:
    unordered_map<entity, int> e2block = { };
    for (auto& e : es) {
        e2block[e] = floor(uniform(random_variable) * et2n_blocks[_e2et[E(e)]]);
    }
    for (auto& [lt, pw] : lt2initial_prob_within) { // TODO: add lt2initial_prob_between
        if (pw > 0) {
            if (!quiet) cout << "  using a block model for \"" << lt << "\"" << endl;
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
    unordered_map<entity, vector<double>> e2coords = { };
    for (auto& [et, dim] : et2dim) {
        for (auto& e : et2es.at(et)) {
            e2coords[e] = vector<double>(dim);
            for (int d=0; d<dim; d++) {
                auto coord = uniform(random_variable);
                e2coords.at(e)[d] = coord;
            }
        }
    }
    for (auto& [lt, ex] : lt2spatial_decay) {
        if (!quiet) cout << "  using a random geometric model for \"" << lt << "\"" << endl;
        assert (lt.rat13 != RT_ID);
        auto et1 = lt.et1, et3 = lt.et3; auto rat13 = lt.rat13; auto dim = et2dim.at(et1);
        assert (et2dim.at(et3) == dim);
        auto es3 = et2es.at(et3);
        for (auto& e1 : et2es.at(et1)) {
            for (auto& e3 : es3) {
                if ((e3 != e1) && ((e3 > e1) || (rat2inv.at(rat13) != rat13))) {
                    double dist2 = 0.0;
                    for (int d=0; d<dim; d++) {
                        dist2 += pow(e2coords.at(e1)[d] - e2coords.at(e3)[d], 2.0);
                    }
                    probability p = exp(-ex * sqrt(dist2));
                    do_random_link(p, e1, rat13, e3);
                }
            }
        }
    }
    if (debug) verify_angle_consistency();
    if (!quiet) cout << "  ...done."<< endl;
}

void init ()
{
    cout << "INITIALIZING..." << endl;
    cout << " MAX_N_INFLT=" << MAX_N_INFLT << ", MAX_N_E=" << MAX_N_E << endl;
    init_randomness();
    init_data();
    init_entities();
    init_relationship_or_action_types();
    init_summary_events();
    init_links();
    init_gexf();
    do_graphviz_diagrams();
    if (debug) {
        dump_data();
        verify_data_consistency();
    }
    cout << "...INITIALIZATION FINISHED." << endl << endl;
}


