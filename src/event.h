// make sure this file is only included once:
#ifndef INC_EVENT_H
#define INC_EVENT_H

/** Performance-critical inline functions for handling of events.
 *
 *  \file
 *
 *  See \ref data_model.h for how a \ref tricl::event relates to other tricl datatypes.
 */

#include <assert.h>

#include "data_model.h"
#include "probability.h"
#include "io.h"
#include "debugging.h"

inline bool event_is_summary(const event& ev) {
    // summary events are encoded by using the negative entity type id as "entity id":
    bool res = ((ev.e1 < 0) || (ev.e3 < 0));
    if (res) {
        assert (ev.e1 < 0);
        assert (ev.e3 < 0);
        assert (ev.ec == EC_EST);
    }
    return res;
}

inline entity_type summary_et1(const event& ev) {
    return (entity_type) -ev.e1;
}
inline entity_type summary_et3(const event& ev) {
    return (entity_type) -ev.e3;
}

/** Return whether a future instance of the event is scheduled.
 */
inline bool event_is_scheduled (
        event& ev, ///< [in] the event to check, passed by reference for performance
        event_data* evd_ ///< [in] the corresponding variable event data, passed by pointer for performance
        )
{
    assert(evd_ == &ev2data.at(ev));
    return (evd_->t > -INFINITY);
}

inline void _schedule_event (event& ev, event_data* evd_, double left_tail, double right_tail)
{
    assert(evd_ == &ev2data.at(ev));
    timepoint t;
    auto ar = evd_->attempt_rate;
    if (ar < 0.0) throw "negative attempt rate";
    assert (ar >= 0.0);
    if (event_is_summary(ev)) { // summary event: use only attempt rate (success will be tested in pop_next_event):
        t = current_t + exponential(random_variable) / (ar * ev2max_success_probability[ev]);
        if (verbose) cout << "         (re)scheduling " << ev << ": summary event, attempt rate " << ar << " → attempt at t=" << t << ", test success then" << endl;
    } else { // particular event: use effective rate:
        auto spu = evd_->success_probunits;
        if (spu == -INFINITY) {
            t = INFINITY;
            if (debug) cout << "         (re)scheduling " << ev << ": zero success probability → t=" << t << endl;
        } else if (ar < INFINITY) {
            auto er = effective_rate(ar, spu, left_tail, right_tail);
            t = current_t + exponential(random_variable) / er;
            if (verbose) {
                if (t==INFINITY) {
                    if (debug) cout << "         (re)scheduling " << ev << ": zero effective rate → t=" << t << endl;
                }
                else if (verbose) cout << "         (re)scheduling " << ev << ": ar " << ar << ", spu " << spu << " → eff. rate " << er << " → next at t=" << t << endl;
            }
        } else {
            // event should happen "right away". to make sure all those events
            // occur in random order, we formally schedule them at some
            // random "past" time instead:
            t = current_t - uniform(random_variable);
            if (verbose) cout << "inf. attempt rate, success probability > 0 → \"immediate\" t=" << t << endl;
        }
    }
    if (t == INFINITY) {
        // replace INFINITY by some unique finite but non-reached time point:
        t = max_t * (1 + uniform(random_variable));
    }
    t2ev[t] = ev;
    evd_->t = t;
}

inline void schedule_event (event& ev, event_data* evd_, double left_tail, double right_tail)
{
    assert(evd_ == &ev2data.at(ev));
    if (event_is_scheduled(ev, evd_)) throw "event already scheduled";
    assert(!event_is_scheduled(ev, evd_));
    _schedule_event(ev, evd_, left_tail, right_tail);
    if (debug) verify_data_consistency();
}

inline void reschedule_event (event& ev, event_data* evd_, double left_tail, double right_tail)
{
    assert(evd_ == &ev2data.at(ev));
    assert(event_is_scheduled(ev, evd_));
    t2ev.erase(evd_->t);
    _schedule_event(ev, evd_, left_tail, right_tail);
    if (debug) verify_data_consistency();
}

void add_event (event& ev);

void remove_event (event& ev, event_data* evd_);

void conditionally_remove_event(event& ev);

void update_adjacent_events (event& ev);

void add_reverse_event (event& old_ev);

void perform_event (event& ev);

bool pop_next_event ();

#endif
