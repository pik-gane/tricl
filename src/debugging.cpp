/** Some stuff only needed for debugging
 *
 *  \file
 */

#include "assert.h"
#include <iostream>

#include "global_variables.h"
#include "angle.h"
#include "io.h"

/** Find the no. of angles influencing an event from scratch
 *  in order to compare it with the stored no.
 */
int compute_n_angles (event_type evt, entity e1, entity e3, bool print) {
    if (print) cout << evt << endl;
    auto outs1 = e2outs[e1];
    auto ins3 = e2ins[e3];
    int na = 0;
    angle_vec as = leg_intersection(e1, outs1, ins3, e3);
    for (auto a_it = as.begin(); a_it < as.end(); a_it++) {
        influence_type inflt = { .evt = evt, .at = { .rat12 = a_it->rat12, .et2 = e2et[a_it->e2], .rat23 = a_it->rat23 } };
        auto dar = _inflt2attempt_rate[INFLT(inflt)];
        auto dsl = _inflt2delta_probunits[INFLT(inflt)];
        if (print) cout << " " << rat2label[a_it->rat12] << " " << e2label[a_it->e2] << " " << rat2label[a_it->rat23] << ", " << INFLT(inflt) << " " << dar << " " << dsl << endl;
        if (COUNT_ALL_ANGLES || (dar != 0.0) || (dsl != 0.0)) { // angle can influence event
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
        assert(n == evd.n_angles);
    }
    // angles -> ec2data: TODO
}

/** Verify that other data is consistent.
 */
void verify_data_consistency () {
    // e2outs:
    for (auto& [e1, outs1] : e2outs) {
        for (auto& [rat13, e3] : outs1) assert (e2ins.at(e3).count({e1, rat13}) == 1);
    }
    // e2ins:
    for (auto& [e3, ins3] : e2ins) {
        for (auto& [e1, rat13] : ins3) assert (e2outs.at(e1).count({rat13, e3}) == 1);
    }
    // ev2data:
    for (auto& [ev, evd] : ev2data) {
        assert (evd.n_angles >= 0);
        assert (evd.attempt_rate >= 0.0);
        assert (evd.success_probunits > -INFINITY);
        if (!(evd.t > -INFINITY)) dump_data();
        assert (evd.t > -INFINITY);
        if (!(((evd.t < INFINITY) && (t2ev.count(evd.t) == 1))
                || ((evd.t == INFINITY) && (t2ev.count(evd.t) > 0))))
            cout << ev << evd << " " << (evd.t == INFINITY) << " " << t2ev.count(evd.t) << endl;
        assert (((evd.t < INFINITY) && (t2ev.count(evd.t) == 1))
                || ((evd.t == INFINITY) && (t2ev.count(evd.t) > 0)));
    }
    // t2be:
    for (auto& [t, ev] : t2ev) {
        assert (t > -INFINITY);
        assert (ev2data.count(ev) == 1);
    }
}
