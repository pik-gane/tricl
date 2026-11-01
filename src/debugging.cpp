/** Some stuff only needed for debugging
 *
 *  \file
 */

#include <iostream>
#include <cstdlib>

/// Consistency check that stays active in optimized builds (unlike assert), used by the --debug mode.
#define CHECK(cond) if (!(cond)) { cerr << "CONSISTENCY CHECK FAILED: " #cond " (" << __FILE__ << ":" << __LINE__ << ")" << endl; abort(); }

#include "global_variables.h"
#include "angle.h"
#include "io.h"
#include "legs.h"
#include "actions.h"

/** Compute total finite event rate from scratch
 *  in order to compare it with the stored one
 */
rate compute_total_finite_er ()
{
    rate ter = 0;
    rate A; probunits B;
    // go through all scheduled events (summary events contribute the rate of all pairs of their types):
    for (auto& [ev, evd] : ev2data)
    {
        auto er = evd.effective_rate;
        if ((er < INFINITY) && !evd.tv) ter += er;  // (time-varying rates are integrated separately, see actions.h)
        if ((ev.ec != EC_TERM) && !event_is_summary(ev))  // event is also covered by summary event
        {
            // subtract the share of this pair that was added when processing the summary event:
            int evt_id = evt_id_of(ev.ec, e2et[ev.e1], ev.rat13, e2et[ev.e3]);
            if ((evt_id >= 0) && evtid2summary[evt_id]) ter -= covered_pair_rate(evt_id, ev.e1, ev.e3, A, B);
        }
    }
    // go through all existing links (their pairs cannot be established):
    for (entity e1 = 1; e1 <= max_e; e1++)
    {
        for (auto& l : e2outs[e1])
        {
            if (l.rat_out == RT_ID) continue;
            int evt_id = evt_id_of(EC_EST, e2et[e1], l.rat_out, e2et[l.e_target]);
            if ((evt_id >= 0) && evtid2summary[evt_id]) ter -= covered_pair_rate(evt_id, e1, l.e_target, A, B);
        }
    }
    // equal entities are not linked:
    for (int evt_id = 0; evt_id < n_evt_ids; evt_id++)
    {
        auto& evt = evtid2evt[evt_id];
        if (!evtid2summary[evt_id] || (evt.et1 != evt.et3)) continue;
        if (!evtid_is_leg_summary[evt_id]) ter -= et2n[evt.et1] * evtid2summary_single_er[evt_id];
        else for (auto& e : et2es[evt.et1]) ter -= covered_pair_rate(evt_id, e, e, A, B);
    }
    return ter;
}

/** Verify that the stored data of every individually scheduled event equals a recomputation from the network,
 *  and that the leg data is consistent.
 */
void verify_rates (bool check_total)
{
    for (auto& [ev, evd] : ev2data)
    {
        if (event_is_summary(ev)) continue;
        int evt_id = evt_id_of(ev.ec, e2et[ev.e1], ev.rat13, e2et[ev.e3]);
        event_data tmp = {};
        tv_data tvd = {};
        bool has_decaying = false;
        compute_event_data(evt_id, ev.e1, ev.rat13, ev.e3, tmp, &tvd, &has_decaying);
        CHECK(has_decaying == (ev2tv.count(ev) > 0));
        if (has_decaying) {
            // the stored weights of decaying influences, decayed to the current time, must equal the recomputed ones:
            const auto& tv = ev2tv.at(ev);
            const auto& dec = evtid2infl_decay[evt_id];
            for (size_t j = 0; j < dec.size(); j++) if (dec[j] > 0) {
                double w = tv.W[j] * exp(- dec[j] * (current_t - tv.t_ref));
                if (!(fabs(w - tvd.W[j]) <= 1e-9 * max(1.0, tvd.W[j]))) cerr << "event " << ev << ": weight of influence " << j << " is " << w << ", recomputed " << tvd.W[j] << endl;
                CHECK(fabs(w - tvd.W[j]) <= 1e-9 * max(1.0, tvd.W[j]));
            }
        }
        for (size_t j = 0; j < MAX_INFL_PER_EVT; j++) if (tmp.n_infl[j] != evd.n_infl[j]) {
            cerr << "event " << ev << ": influence " << j << " counted " << evd.n_infl[j] << " times, recomputed " << tmp.n_infl[j] << endl;
            CHECK(tmp.n_infl[j] == evd.n_infl[j]);
        }
        CHECK(fabs(tmp.attempt_rate - evd.attempt_rate) <= 1e-9 * max(1.0, fabs(tmp.attempt_rate)));
        CHECK(fabs(tmp.success_probunits - evd.success_probunits) <= 1e-9 * max(1.0, fabs(tmp.success_probunits)));
        CHECK((tmp.n_inf_attempt == evd.n_inf_attempt) && (tmp.n_pos_inf_probunits == evd.n_pos_inf_probunits) && (tmp.n_neg_inf_probunits == evd.n_neg_inf_probunits));
        if (track_legs) CHECK((e2src_evs[ev.e1].count(ev) == 1) && (e2tgt_evs[ev.e3].count(ev) == 1));
    }
    if (check_total) {
        // the incrementally maintained total rate (used by the log-likelihood):
        rate exact = compute_total_finite_er();
        if (!(fabs(exact - total_finite_effective_rate) <= 1e-8 * max(1.0, fabs(exact))))
            cerr << "total rate " << total_finite_effective_rate << " differs from recomputed " << exact << endl;
        CHECK(fabs(exact - total_finite_effective_rate) <= 1e-8 * max(1.0, fabs(exact)));
    }
    if (!track_legs) return;
    // leg degrees:
    for (entity e = 1; e <= max_e; e++)
    {
        vector<int> deg(legkeys.size(), 0);
        for (auto& l : e2outs[e]) if (l.rat_out != RT_ID) { int k = slot2legkey[at_slot(l.rat_out, e2et[l.e_target], NO_RAT)]; if (k >= 0) deg[k]++; }
        for (auto& l : e2ins[e]) if (l.rat_in != RT_ID) { int k = slot2legkey[at_slot(NO_RAT, e2et[l.e_source], l.rat_in)]; if (k >= 0) deg[k]++; }
        for (size_t k = 0; k < legkeys.size(); k++) CHECK(deg[k] == legdeg(e, k));
        size_t n = 0;
        for (auto& ev : e2src_evs[e]) { CHECK(ev2data.count(ev) == 1); n++; }
        for (auto& ev : e2tgt_evs[e]) CHECK(ev2data.count(ev) == 1);
        (void) n;
    }
    // histograms and totals of summary events with legs:
    for (int evt_id = 0; evt_id < n_evt_ids; evt_id++)
    {
        if (!evtid_is_leg_summary[evt_id]) continue;
        auto& evt = evtid2evt[evt_id];
        map<vector<int>, long> H1, H3;
        for (auto& e : et2es[evt.et1]) H1[outdegs(evt_id, e)]++;
        for (auto& e : et2es[evt.et3]) H3[indegs(evt_id, e)]++;
        CHECK(H1 == evtid2lsd[evt_id].H1);
        CHECK(H3 == evtid2lsd[evt_id].H3);
        double T = recompute_leg_summary_T(evt_id);
        CHECK(fabs(T - evtid2lsd[evt_id].T) <= 1e-9 * max(1.0, T));
    }
}

/** Find the no. of angles influencing an event from scratch
 *  in order to compare it with the stored no.
 */
int compute_n_angles (event_type evt, entity e1, entity e3, bool print) {
    if (print) cout << evt << endl;
    const auto& outs1 = e2outs[e1];
    const auto& ins3 = e2ins[e3];
    int na = 0;
    angle_vec as;
    get_angles(e1, outs1, ins3, e3, as);
    int evt_id = evt_id_of(evt.ec, evt.et1, evt.rat13, evt.et3);
    for (auto a_it = as.begin(); a_it < as.end(); a_it++) {
        int j = (evt_id >= 0) ? infl_index_of(evt_id, a_it->rat12, e2et[a_it->e2], a_it->rat23) : -1;
        if (print) cout << " " << rat2label[a_it->rat12] << " " << e2label[a_it->e2] << " " << rat2label[a_it->rat23] << ", influence index " << j << endl;
        if (j >= 0) { // angle can influence event
            na++;
        }
    }
    return na;
}

/** Verify that data about angles is consistent.
 */
void verify_angle_consistency (bool at_end_of_event) {
    verify_rates(at_end_of_event);
    // ev2data -> n_angles:
    for (auto& [ev, evd] : ev2data) {
        if (event_is_summary(ev)) continue;
        auto e1 = ev.e1, e3=ev.e3;
        auto et1 = e2et[e1], et3 = e2et[e3];
        auto n = compute_n_angles({ ev.ec, et1, ev.rat13, et3 }, e1, e3, false);
        if (n != evd.n_angles) {
           cout << "failed at " << ev << " " << evd << " " << n << endl;
           compute_n_angles({ ev.ec, et1, ev.rat13, et3 }, e1, e3, true);
//           dump_data();
        }
        CHECK(n == evd.n_angles);
    }
    // angles -> ec2data: TODO
}

/** Verify that other data is consistent.
 */
void verify_data_consistency () {
    // e2outs:
    for (entity e1 = 1; e1 <= max_e; e1++) {
        for (auto& l : e2outs[e1]) {
            auto rat13 = l.rat_out;
            auto e3 = l.e_target;
            CHECK(e2ins.at(e3).count({e1, rat13}) == 1);
        }
    }
    // e2ins:
    for (entity e3 = 1; e3 <= max_e; e3++) {
        for (auto& l : e2ins[e3]) {
            auto e1 = l.e_source;
            auto rat13 = l.rat_in;
            CHECK(e2outs.at(e1).count({rat13, e3}) == 1);
        }
    }
    // ev2data:
    for (auto& [ev, evd] : ev2data) {
        CHECK(evd.n_angles >= 0);
        CHECK(evd.attempt_rate >= 0.0);
        CHECK(evd.attempt_rate < INFINITY);
        CHECK(std::isfinite(evd.success_probunits));
        CHECK((evd.n_inf_attempt >= 0) && (evd.n_pos_inf_probunits >= 0) && (evd.n_neg_inf_probunits >= 0));
        if (!(evd.t > -INFINITY)) dump_data();
        CHECK(evd.t > -INFINITY);
        if (scheduling_enabled) {
            // the event is either in the rate tree (finite rate) or in the list of immediate events (infinite rate):
            CHECK((evd.slot >= 0) != (evd.imm >= 0));
            if (evd.slot >= 0) {
                CHECK(schedule.is_used(evd.slot) && (schedule.at(evd.slot) == ev));
                CHECK(evd.effective_rate < INFINITY);
                if (!event_is_summary(ev)) CHECK(schedule.weight(evd.slot) == evd.effective_rate);
            } else {
                CHECK((evd.imm < (int) immediate_events.size()) && (immediate_events[evd.imm] == ev));
                CHECK(evd.effective_rate == INFINITY);
            }
        } else {
            CHECK((evd.slot < 0) && (evd.imm < 0));
        }
    }
    // the schedule:
    if (scheduling_enabled) {
        CHECK(schedule.size() + (int) immediate_events.size() == (int) ev2data.size());
        double exact = schedule.exact_total();
        CHECK(fabs(schedule.total() - exact) <= 1e-9 * max(1.0, exact));
        for (auto& ev : immediate_events) CHECK(ev2data.count(ev) == 1);
    } else {
        CHECK((schedule.size() == 0) && immediate_events.empty());
    }
}
