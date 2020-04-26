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
    rate ar = evd_->attempt_rate;
    if (ar < 0.0) throw "negative attempt rate";
    auto spu = evd_->success_probunits;
    timepoint t;
    if (event_is_summary(ev))  // summary event:
    {
        // use a common upper bound to the actual effective rate for scheduling (actual success will then later be tested in pop_next_event):
        t = current_t + exponential(random_variable) / (ar * summary_ev2max_success_probability[ev]);
        if (verbose) cout << "         (re)scheduling " << ev << ": summary event, attempt rate " << ar << " → attempt at t=" << t << ", test success then" << endl;
        // compute base effective rate using base success probability units:
        rate er = evd_->effective_rate = effective_rate(ar, spu, left_tail, right_tail);
        assert (er < INFINITY);
        // register it in total:
        add_effective_rate(er);
        // NOTE: in update_adjacent_events, particular events must get update total effective rate properly as the difference to this!
    }
    else  // particular event:
    {
        // use effective rate for scheduling:
        if (spu == -INFINITY)
        {
            evd_->effective_rate = 0;
            t = INFINITY;
            if (debug) cout << "         (re)scheduling " << ev << ": zero success probability → t=" << t << endl;
        }
        else if (ar < INFINITY)
        {
            // compute effective rate:
            rate er = evd_->effective_rate = effective_rate(ar, spu, left_tail, right_tail);
            assert (er < INFINITY);
            // register it in total:
            add_effective_rate(er);

            // draw time interval after which it would happen if nothing changes in between:
            timepoint dt = exponential(random_variable) / er;
            // add it to current time to get occurence time:
            t = current_t + dt;

            if (verbose) {
                if (t==INFINITY) {
                    if (debug) cout << "         (re)scheduling " << ev << ": zero effective rate → t=" << t << endl;
                }
                else if (verbose) cout << "         (re)scheduling " << ev << ": ar " << ar << ", spu " << spu << " → eff. rate " << er << " → next at t=" << t << endl;
            }
        }
        else  // event should happen "right away"
        {
            // to make sure all those events occur in random order,
            // we formally schedule them at some random "past" time instead:
            rate er = evd_->effective_rate = INFINITY;
            // register it in total:
            add_effective_rate(er);
            t = current_t - abs(1 + current_t) * uniform(random_variable);
            if (verbose) cout << "         (re)scheduling " << ev << ": ar inf, spu > 0 → eff. rate inf → next \"immediately\" at t=" << t << endl;
        }
    }
    if (t == INFINITY)
    {
        // replace INFINITY by some unique finite but non-reached time point:
        t = max_t * (1 + uniform(random_variable));
    }
    // store time:
    evd_->t = t;
    t2ev[t] = ev;
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

    // remove from schedule and total_effective_rate:
    t2ev.erase(evd_->t);
    subtract_effective_rate(evd_->effective_rate, !event_is_summary(ev));

    // schedule anew:
    _schedule_event(ev, evd_, left_tail, right_tail);
    if (debug) verify_data_consistency();
}

void add_event (event& ev);

void remove_event (event& ev, event_data* evd_);

void conditionally_remove_event(event& ev);

void update_adjacent_events (event& ev);

void add_reverse_event (event& old_ev);

void perform_event (event& ev, event_data* evd_);

bool pop_next_event ();

#endif
