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

/** Compute total finite event rate from scratch
 *  in order to compare it with the stored one
 */
rate compute_total_finite_er ()
{
    rate ter = 0;
    // go through all scheduled events:
    for (auto& [ev, evd] : ev2data)
    {
        auto ec = ev.ec;
        auto er = evd.effective_rate;
        if (er < INFINITY) ter += er;
        if ((ec != EC_TERM) && !event_is_summary(ev))  // event is also covered by summary event
        {
            // subtract single er that will be added when processing summary event:
            ter -= summary_single_er_of(ec, e2et[ev.e1], ev.rat13, e2et[ev.e3]);
        }
    }
    // go through all existing links:
    for (entity e1 = 1; e1 <= max_e; e1++)
    {
        for (auto& l : e2outs[e1])
        {
            // subtract single er that was added when processing summary event:
            ter -= summary_single_er_of(EC_EST, e2et[e1], l.rat_out, e2et[l.e_target]);
        }
    }
    // go through all entity types:
    for (auto& [et, n] : et2n)
    {
        for (auto& rat13 : ets2relations[{et, et}])
        {
            // subtract n times single er since equal entities will not be linked:
            ter -= n * summary_single_er_of(EC_EST, et, rat13, et);
        }
    }

    return ter;
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
void verify_angle_consistency () {
    // ev2data -> n_angles:
    for (auto& [ev, evd] : ev2data) {
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
