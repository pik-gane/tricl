/** Do initialization tasks before simulation, after configuration.
 *
 *  \file
 */

#include <sstream>

#include "global_variables.h"
#include "probability.h"
#include "entity.h"
#include "link.h"
#include "event.h"
#include "graphviz.h"
#include "gexf.h"
#include "init.h"
#include "io.h"
#include "debugging.h"

// parameters:
int n_rats = 0; // total no. of rats
unordered_set<event_type> possible_evts = {};
bool any_leg_influences = false;
entity_type e2et[MAX_N_E];

// dense type indexing (see init_types()):
int n_et_slots = 0, n_rat_slots = 0, n_at_slots = 0, n_evt_ids = 0;
vector<int> evt_slot2id;
vector<event_type> evtid2evt;
vector<rate> evtid2base_attempt_rate;
vector<probunits> evtid2base_probunits;
vector<double> evtid2left_tail, evtid2right_tail, evtid2scale;
vector<rate> evtid2summary_single_er;
vector<probability> evtid2summary_max_success_probability;
vector<rate> inflt_attempt_rate;
vector<probunits> inflt_delta_probunits;
vector<vector<relationship_or_action_type>> ets2rats;
vector<vector<int>> evtid2infl_slots;
vector<signed char> evtid_at2infl;

// model parameters and gradients:
vector<model_param> params;
vector<int> evtid2base_attempt_param, evtid2base_probunits_param;
vector<vector<int>> evtid2infl_attempt_param, evtid2infl_probunits_param;
vector<double> grad_event_terms, grad_rate, grad_exposure;

// derived constants:
unordered_set<entity> es;
unordered_map<entity_type_pair, unordered_set<relationship_or_action_type>> ets2relations;  // possible relations

// variable data:

timepoint current_t = 0, last_dt = 0;
long int n_events = 0;
event current_ev = {};
event_data* current_evd_ = &sure_evd;

// network state:
unordered_map<entity_type, vector<entity>> et2es = {};  // kept to equal inverse of e2et
vector<outleg_set> e2outs = {};
vector<inleg_set> e2ins = {};
unordered_map<link_type, long int> lt2n = {};
long int n_links = 0, n_angles = 0;

// event data:
unordered_map<event, event_data> ev2data = {};
map<timepoint, event> t2ev = {};  // kept to equal inverse of ev2data.t

// log-likelihood:
double cumulative_logl = 0;
rate total_finite_effective_rate = 0;
int n_infinite_effective_rates = 0;


/** \returns the label of an event type as printed by operator<<. */
static string evt_to_string (const event_type& evt)
{
    std::ostringstream os;
    os << evt;
    return os.str();
}

/** Look up a value in a map of per-event-type parameters, returning a default if missing. */
template <class M, class V>
static V evt_param (const M& m, const event_type& evt, V deflt)
{
    auto it = m.find(evt);
    return (it == m.end()) ? deflt : it->second;
}

/** Set up the dense type indexing and the parameter tables used in the simulation's hot loops
 *  (see global_variables.h), from the maps filled by read_config().
 *
 *  Must be called after init_relationship_or_action_types() and before init_events().
 */
void init_types ()
{
    // slots:
    n_et_slots = 1;
    for (auto& [et, l] : et2label) n_et_slots = max(n_et_slots, (int) et + 1);
    n_rat_slots = 1;
    for (auto& [rat, l] : rat2label) n_rat_slots = max(n_rat_slots, (int) rat + 1);
    n_at_slots = n_rat_slots * n_et_slots * n_rat_slots;

    // possible event types get dense ids:
    for (auto& [evt, ar] : evt2base_attempt_rate) {
        if (ar > 0.0) possible_evts.insert(evt);
    }
    for (auto& [inflt, ar] : inflt2attempt_rate) {
        if (ar > 0.0) possible_evts.insert(inflt.evt);
    }
    evt_slot2id.assign(3 * n_et_slots * n_rat_slots * n_et_slots, -1);
    evtid2evt.clear();
    n_evt_ids = 0;
    for (auto& evt : possible_evts) {
        evt_slot2id[evt_slot(evt.ec, evt.et1, evt.rat13, evt.et3)] = n_evt_ids++;
        evtid2evt.push_back(evt);
    }

    // per-event-type parameters:
    evtid2base_attempt_rate.assign(n_evt_ids, 0.0);
    evtid2base_probunits.assign(n_evt_ids, 0.0);
    evtid2left_tail.assign(n_evt_ids, 1.0);
    evtid2right_tail.assign(n_evt_ids, 1.0);
    evtid2scale.assign(n_evt_ids, 0.0);
    evtid2summary_single_er.assign(n_evt_ids, 0.0);
    evtid2summary_max_success_probability.assign(n_evt_ids, 0.0);
    for (int id = 0; id < n_evt_ids; id++) {
        auto& evt = evtid2evt[id];
        rate ar1 = evtid2base_attempt_rate[id] = evt_param(evt2base_attempt_rate, evt, (rate) 0.0);
        probunits spu0 = evtid2base_probunits[id] = evt_param(evt2base_probunits, evt, (probunits) 0.0);
        double left_tail = evtid2left_tail[id] = evt_param(evt2left_tail, evt, 1.0),
               right_tail = evtid2right_tail[id] = evt_param(evt2right_tail, evt, 1.0),
               scale = evtid2scale[id] = tail2scale(left_tail) + tail2scale(right_tail);
        if ((evt.ec == EC_EST) && (ar1 > 0.0)) {
            // this event type has a summary event for its spontaneous occurrence:
            evtid2summary_single_er[id] = effective_rate(ar1, spu0, left_tail, right_tail, scale);
            // compile maximal success units. if no influences can increase the success units,
            // this equals the base_probunits, otherwise it is infinite:
            probunits max_spu = spu0;
            for (auto& [inflt, pu] : inflt2delta_probunits) {
                if ((inflt.evt == evt) && (pu > 0.0)) {
                    max_spu = INFINITY;
                    break;
                }
            }
            evtid2summary_max_success_probability[id] = probunits2probability(max_spu, left_tail, right_tail, scale);
        }
    }

    // influence tables:
    inflt_attempt_rate.assign((size_t) n_evt_ids * n_at_slots, 0.0);
    inflt_delta_probunits.assign((size_t) n_evt_ids * n_at_slots, 0.0);
    any_leg_influences = false;
    for (auto& [inflt, ar] : inflt2attempt_rate) {
        assert (!( inflt.evt.ec == EC_EST && ( inflt.at.rat12 == NO_RAT || inflt.at.rat23 == NO_RAT ) ));
        int id = evt_id_of(inflt.evt.ec, inflt.evt.et1, inflt.evt.rat13, inflt.evt.et3);
        if (id >= 0) inflt_attempt_rate[inflt_index(id, inflt.at.rat12, inflt.at.et2, inflt.at.rat23)] = ar;
        if (((inflt.at.rat12 == NO_RAT) || (inflt.at.rat23 == NO_RAT)) && (ar != 0.0)) any_leg_influences = true;
    }
    for (auto& [inflt, spu] : inflt2delta_probunits) {
        assert (!( inflt.evt.ec == EC_EST && ( inflt.at.rat12 == NO_RAT || inflt.at.rat23 == NO_RAT ) ));
        int id = evt_id_of(inflt.evt.ec, inflt.evt.et1, inflt.evt.rat13, inflt.evt.et3);
        if (id >= 0) inflt_delta_probunits[inflt_index(id, inflt.at.rat12, inflt.at.et2, inflt.at.rat23)] = spu;
        if (((inflt.at.rat12 == NO_RAT) || (inflt.at.rat23 == NO_RAT)) && (spu != 0.0)) any_leg_influences = true;
    }

    // possible relationship or action types by entity type pair
    // (in the same order as in ets2relations, since this order determines the order of random draws):
    ets2rats.assign((size_t) n_et_slots * n_et_slots, {});
    for (auto& [ets, rats] : ets2relations) {
        for (auto& rat : rats) ets2rats[(int) ets.et1 * n_et_slots + (int) ets.et3].push_back(rat);
    }

    // influences (angle types with nonzero effect) by event type, and the model parameters:
    evtid2infl_slots.assign(n_evt_ids, {});
    evtid_at2infl.assign((size_t) n_evt_ids * n_at_slots, -1);
    params.clear();
    evtid2base_attempt_param.assign(n_evt_ids, -1);
    evtid2base_probunits_param.assign(n_evt_ids, -1);
    evtid2infl_attempt_param.assign(n_evt_ids, {});
    evtid2infl_probunits_param.assign(n_evt_ids, {});
    auto register_infl = [&](int id, const angle_type& at) {
        int slot = at_slot(at.rat12, at.et2, at.rat23);
        if (evtid_at2infl[(size_t) id * n_at_slots + slot] >= 0) return;  // already registered
        if ((inflt_attempt_rate[(size_t) id * n_at_slots + slot] == 0.0) && (inflt_delta_probunits[(size_t) id * n_at_slots + slot] == 0.0)) return;  // no effect
        int j = evtid2infl_slots[id].size();
        if (j >= MAX_INFL_PER_EVT) throw "too many angle types influence event type \"" + evt_to_string(evtid2evt[id]) + "\" (recompile with larger MAX_INFL_PER_EVT)";
        evtid2infl_slots[id].push_back(slot);
        evtid_at2infl[(size_t) id * n_at_slots + slot] = j;
        evtid2infl_attempt_param[id].push_back(-1);
        evtid2infl_probunits_param[id].push_back(-1);
    };
    // (iterate in a deterministic order: by event type id, then by the config maps)
    for (int id = 0; id < n_evt_ids; id++) {
        auto& evt = evtid2evt[id];
        for (auto& [inflt, ar] : inflt2attempt_rate) if (inflt.evt == evt) register_infl(id, inflt.at);
        for (auto& [inflt, spu] : inflt2delta_probunits) if (inflt.evt == evt) register_infl(id, inflt.at);
        // parameters (only those specified in the config file; influence parameters only for influences with nonzero effect):
        string evt_label = evt_to_string(evt);
        if (evt2base_attempt_rate.count(evt) > 0) {
            evtid2base_attempt_param[id] = params.size();
            params.push_back({ PK_BASE_ATTEMPT, id, -1, evtid2base_attempt_rate[id], evt_label + " | base attempt" });
        }
        if (evt2base_probunits.count(evt) > 0) {
            evtid2base_probunits_param[id] = params.size();
            params.push_back({ PK_BASE_PROBUNITS, id, -1, evtid2base_probunits[id], evt_label + " | base probunits" });
        }
        for (size_t j = 0; j < evtid2infl_slots[id].size(); j++) {
            int slot = evtid2infl_slots[id][j];
            // reconstruct the angle type from the slot:
            angle_type at = { .rat12 = (relationship_or_action_type) (slot / (n_et_slots * n_rat_slots)),
                              .et2 = (entity_type) ((slot / n_rat_slots) % n_et_slots),
                              .rat23 = (relationship_or_action_type) (slot % n_rat_slots) };
            string at_label = rat2label[at.rat12] + " " + et2label[at.et2] + " " + rat2label[at.rat23];
            influence_type inflt = { .evt = evt, .at = at };
            if (inflt2attempt_rate.count(inflt) > 0) {
                evtid2infl_attempt_param[id][j] = params.size();
                params.push_back({ PK_INFL_ATTEMPT, id, (int) j, inflt_attempt_rate[(size_t) id * n_at_slots + slot], evt_label + " | attempt via " + at_label });
            }
            if (inflt2delta_probunits.count(inflt) > 0) {
                evtid2infl_probunits_param[id][j] = params.size();
                params.push_back({ PK_INFL_PROBUNITS, id, (int) j, inflt_delta_probunits[(size_t) id * n_at_slots + slot], evt_label + " | probunits via " + at_label });
            }
        }
    }
    grad_event_terms.assign(params.size(), 0.0);
    grad_rate.assign(params.size(), 0.0);
    grad_exposure.assign(params.size(), 0.0);

    if (verbose) cout << " " << n_evt_ids << " possible event types, " << n_at_slots << " angle type slots, " << params.size() << " model parameters" << endl;
}

/** Prepare all entities.
 */
void init_entities ()
{
    // inspect pre-registered entities:
    auto et2remaining_n = et2n;
    for (auto& [e, l] : e2label) {
        auto et = e2et[e];
        assert(e >= 0);
        if (et2remaining_n[et] > 0) {
            et2remaining_n[et]--;
        } else {
            et2n[et]++;
        }
    }

    // generate remaining entities:
    for (auto& [et, n] : et2remaining_n) {
        if (!silent) cout << " entity type \"" << et2label[et] << "\" has " << et2n[et] << " entities" << endl;
        assert (n >= 0);
        while (n > 0) {
            add_entity(et, "");
            n--;
        }
        assert ((entity) et2es[et].size() == et2n[et]);
    }
    if (max_e >= 1<<E_BITS) throw "too many entities (recompile with larger E_BITS?)";
}

/** Analyse relationship or action types.
 */
void init_relationship_or_action_types ()
{
    // verify symmetry of relationship inversion map:
    assert (rat2inv.at(RT_ID) == RT_ID);
    for (auto& [r, la] : rat2label) {
        if (rat2inv.count(r) == 0) rat2inv[r] = NO_RAT;
        auto inv = rat2inv[r];
        if ((inv != NO_RAT) && (inv != r)) {
            assert (rat2inv.at(inv) == r);
            assert (r_is_action_type.at(r) == r_is_action_type.at(inv));
        }
    }

    // register possible relationship types by entity type pair, and compute probunits:
    for (auto& [evt, pu] : evt2base_probunits) {
        auto rat13 = evt.rat13;
        assert (rat13 != RT_ID);
        entity_type_pair ets = { evt.et1, evt.et3 };
        if (ets2relations.count(ets) == 0) ets2relations[ets] = {};
        ets2relations[ets].insert(rat13);
        lt2n[{evt.et1, rat13, evt.et3}] = 0;
    }

    n_rats = rat2label.size();
}

/** Set up initial schedule of summary events.
 */
void init_events ()
{
    if (verbose) {
        if (!silent) cout << " possible event types with base attempt rates and base success probabilities:" << endl;
        for (int id = 0; id < n_evt_ids; id++) cout << "  " << evtid2evt[id] << ": " << evtid2base_attempt_rate[id] <<
                ", " << probunits2probability(evtid2base_probunits[id], evtid2left_tail[id], evtid2right_tail[id], evtid2scale[id]) << endl;
    }

    // summary events for purely spontaneous establishment without angles:

    if (!quiet) cout << " initial scheduling of summary events..." << endl;
    for (auto& [ets, relations] : ets2relations) {
        auto et1 = ets.et1, et3 = ets.et3;
        for (auto& rat13 : relations) {
            int evt_id = evt_id_of(EC_EST, et1, rat13, et3);
            if (evt_id < 0) continue;
            auto ar1 = evtid2base_attempt_rate[evt_id];
            if (ar1 > 0) {
                event summary_ev = { .ec = EC_EST,
                        .e1 = (entity)-et1, // in spontaneous events, fields e1 and e3 are used to store entity types with negative sign
                        .rat13 = rat13,
                        .e3 = (entity)-et3 };
                if (verbose) cout << "  " << et2label[et1] << " " << rat2label[rat13] << " " << et2label[et3] << endl;
                rate ar_all = ar1 * et2n[et1] * et2n[et3];
                event_data summary_evd = {};
                summary_evd.attempt_rate = ar_all;  // finite, see config validation
                add_probunits_contribution(&summary_evd, evtid2base_probunits[evt_id]);
                summary_evd.t = -INFINITY;
                ev2data[summary_ev] = summary_evd;
                schedule_event(summary_ev, &ev2data[summary_ev], evt_id);
                // (the gradient contributions of the summary event are registered once here since its rate never changes)
                register_summary_share_gradient(evt_id, (double) et2n[et1] * et2n[et3]);
                // adjust effective rate because equal entities won't be linked
                // (each of the n excluded pairs had contributed the single effective rate ar1 * p0):
                if (et1 == et3)
                {
                    subtract_summary_shares(evt_id, et2n[et1]);
                }
            }
        }
    }
    if (!quiet) cout << "  ...done."<< endl;
}

/** Add all initial links and corresponding events.
 */
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
        event ev = { .ec=EC_EST, l.e1, l.rat13, l.e3 };
        conditionally_remove_event(ev);
        // prepair total effective rate since call to perform_event will subtracted INFINITY from it:
        perform_event(ev, &sure_evd);  // this automatically also adds the inverse link, if any.
    }

    // random links:

    // block model:
    unordered_map<entity, int> e2block = { };
    for (auto& e : es) {
        e2block[e] = floor(uniform(initial_state_random_variable) * et2n_blocks[e2et[e]]);
    }
    for (auto& [lt, pw] : lt2initial_prob_within) { // TODO: add lt2initial_prob_between
        if (pw > 0) {
            if (verbose) cout << "  using a block model for \"" << lt << "\"" << endl;
            assert (lt.rat13 != RT_ID);
            auto et1 = lt.et1, et3 = lt.et3; auto rat13 = lt.rat13;
            for (auto& e1 : et2es[et1]) {
                for (auto& e3 : et2es[et3]) {
                    // (for symmetric relationship types between entities of the same type, each pair is considered once)
                    if ((e3 != e1) && ((et1 != et3) || (e3 > e1) || (rat2inv[rat13] != rat13))) {
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
                auto coord = uniform(initial_state_random_variable);
                e2coords.at(e)[d] = coord;
            }
        }
    }
    for (auto& [lt, ex] : lt2spatial_decay) {
        if (verbose) cout << "  using a random geometric model for \"" << lt << "\"" << endl;
        assert (lt.rat13 != RT_ID);
        auto et1 = lt.et1, et3 = lt.et3; auto rat13 = lt.rat13; auto dim = et2dim.at(et1);
        assert (et2dim.at(et3) == dim);
        auto es3 = et2es.at(et3);
        for (auto& e1 : et2es.at(et1)) {
            for (auto& e3 : es3) {
                if ((e3 != e1) && ((et1 != et3) || (e3 > e1) || (rat2inv.at(rat13) != rat13))) {
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
    // reset cumulative loglikelihood (and its gradient) to count only what happens after initial state:
    cumulative_logl = 0;
    std::fill(grad_event_terms.begin(), grad_event_terms.end(), 0.0);
    std::fill(grad_exposure.begin(), grad_exposure.end(), 0.0);

    if (debug) verify_angle_consistency();
    if (!quiet) cout << "  ...done."<< endl;
}

/** Only set up the type indexing and model parameters (for --dump-parameters).
 */
void init_parameters_only ()
{
    init_relationship_or_action_types();
    init_types();
}

/** Perform all initialization tasks.
 */
void init ()
{
    if (!silent) cout << "INITIALIZING..." << endl;
    init_randomness();
    init_entities();
    init_relationship_or_action_types();
    init_types();
    init_events();
    init_links();
    open_events_out();  // only after initial links, so that only simulated events are written
    open_stats_out();
    open_links_out();
    write_entities_out();
    init_gexf();
    do_graphviz_diagrams();
    if (debug) {
        dump_data();
        verify_data_consistency();
    }
    if (!silent) cout << "...INITIALIZATION FINISHED." << endl << endl;
}

