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

// Helpers for accumulating attempt rate and success probunits contributions
// (infinite contributions are counted rather than summed, see event_data):

/** Add a contribution to the attempt rate of an event. */
inline void add_attempt_contribution (event_data* evd_, rate dar)
{
    if (dar == INFINITY) evd_->n_inf_attempt++;
    else evd_->attempt_rate += dar;
}
/** Remove a contribution from the attempt rate of an event. */
inline void remove_attempt_contribution (event_data* evd_, rate dar)
{
    if (dar == INFINITY) { evd_->n_inf_attempt--; assert (evd_->n_inf_attempt >= 0); }
    else evd_->attempt_rate = max(0.0, evd_->attempt_rate - dar);
}
/** Add a contribution to the success probunits of an event. */
inline void add_probunits_contribution (event_data* evd_, probunits dspu)
{
    if (dspu == INFINITY) evd_->n_pos_inf_probunits++;
    else if (dspu == -INFINITY) evd_->n_neg_inf_probunits++;
    else evd_->success_probunits += dspu;
}
/** Remove a contribution from the success probunits of an event. */
inline void remove_probunits_contribution (event_data* evd_, probunits dspu)
{
    if (dspu == INFINITY) { evd_->n_pos_inf_probunits--; assert (evd_->n_pos_inf_probunits >= 0); }
    else if (dspu == -INFINITY) { evd_->n_neg_inf_probunits--; assert (evd_->n_neg_inf_probunits >= 0); }
    else evd_->success_probunits -= dspu;
}
/** \returns the total attempt rate of an event, taking infinite contributions into account. */
inline rate total_attempt_rate (const event_data* evd_)
{
    return (evd_->n_inf_attempt > 0) ? INFINITY : evd_->attempt_rate;
}
/** \returns the total success probunits of an event, taking infinite contributions into account.
 *  A -inf contribution (impossible) dominates a +inf contribution (certain).
 */
inline probunits total_success_probunits (const event_data* evd_)
{
    return (evd_->n_neg_inf_probunits > 0) ? -INFINITY : (evd_->n_pos_inf_probunits > 0) ? INFINITY : evd_->success_probunits;
}

// Gradient bookkeeping (only active if compute_gradient is set).
//
// The log-likelihood is sum_i log(rate of event i) - integral of the total rate over time.
// For an event with attempt rate A = a0 + sum_j n_j a_j and success probability sigma(B), B = b0 + sum_j n_j b_j,
// where n_j is the number of adjacent angles of influence j, the rate is A * sigma(B), so that
//   d rate / d a0 = sigma,  d rate / d a_j = n_j sigma,  d rate / d b0 = A sigma',  d rate / d b_j = n_j A sigma',
// and the gradient of log(rate) is the same divided by the rate.
// grad_rate holds the gradient of the total rate and is updated whenever an event's rate is registered or unregistered;
// its integral over time (grad_exposure) is accumulated in advance_time().

/** Add (sign = +1) or remove (sign = -1) the gradient contributions of a particular event's current rate to \ref grad_rate.
 */
inline void register_event_rate_gradient (const event_data* evd_, int evt_id, double sign)
{
    if (!compute_gradient) return;
    rate A = total_attempt_rate(evd_);
    if (!(A < INFINITY) || (A == 0)) return;  // immediate or impossible events have no (finite) gradient
    probunits B = total_success_probunits(evd_);
    if (B == -INFINITY) return;
    double sigma = probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]),
           dsigma = probunits2probability_derivative(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
    int k = evtid2base_attempt_param[evt_id];
    if (k >= 0) grad_rate[k] += sign * sigma;
    const auto& ia = evtid2infl_attempt_param[evt_id];
    for (size_t j = 0; j < ia.size(); j++) {
        if ((ia[j] >= 0) && (evd_->n_infl[j] > 0)) grad_rate[ia[j]] += sign * evd_->n_infl[j] * sigma;
    }
    if (dsigma != 0.0) {
        k = evtid2base_probunits_param[evt_id];
        if (k >= 0) grad_rate[k] += sign * A * dsigma;
        const auto& ib = evtid2infl_probunits_param[evt_id];
        for (size_t j = 0; j < ib.size(); j++) {
            if ((ib[j] >= 0) && (evd_->n_infl[j] > 0)) grad_rate[ib[j]] += sign * evd_->n_infl[j] * A * dsigma;
        }
    }
}

/** Add the gradient contributions of a number of pairs covered by the summary event of an event type to \ref grad_rate
 *  (count may be negative to remove pairs).
 */
inline void register_summary_share_gradient (int evt_id, double count)
{
    if (!compute_gradient || (count == 0.0)) return;
    rate a = evtid2base_attempt_rate[evt_id];
    probunits b = evtid2base_probunits[evt_id];
    if (b == -INFINITY) return;
    int k = evtid2base_attempt_param[evt_id];
    if (k >= 0) grad_rate[k] += count * probunits2probability(b, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
    k = evtid2base_probunits_param[evt_id];
    if (k >= 0) grad_rate[k] += count * a * probunits2probability_derivative(b, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
}

/** Add the gradient of log(rate) of a particular event that happens now to \ref grad_event_terms.
 */
inline void add_event_log_gradient (const event_data* evd_, int evt_id)
{
    if (!compute_gradient) return;
    rate A = total_attempt_rate(evd_);
    if (!(A < INFINITY) || (A == 0)) return;
    probunits B = total_success_probunits(evd_);
    if (B == -INFINITY) return;
    double sigma = probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]),
           dsigma = probunits2probability_derivative(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
    int k = evtid2base_attempt_param[evt_id];
    if (k >= 0) grad_event_terms[k] += 1 / A;
    const auto& ia = evtid2infl_attempt_param[evt_id];
    for (size_t j = 0; j < ia.size(); j++) {
        if ((ia[j] >= 0) && (evd_->n_infl[j] > 0)) grad_event_terms[ia[j]] += evd_->n_infl[j] / A;
    }
    if ((dsigma != 0.0) && (sigma > 0.0)) {
        k = evtid2base_probunits_param[evt_id];
        if (k >= 0) grad_event_terms[k] += dsigma / sigma;
        const auto& ib = evtid2infl_probunits_param[evt_id];
        for (size_t j = 0; j < ib.size(); j++) {
            if ((ib[j] >= 0) && (evd_->n_infl[j] > 0)) grad_event_terms[ib[j]] += evd_->n_infl[j] * dsigma / sigma;
        }
    }
}

/** Add the gradient of log(rate) of an event that happens now via a summary event to \ref grad_event_terms.
 */
inline void add_summary_event_log_gradient (int evt_id)
{
    if (!compute_gradient) return;
    rate a = evtid2base_attempt_rate[evt_id];
    probunits b = evtid2base_probunits[evt_id];
    int k = evtid2base_attempt_param[evt_id];
    if ((k >= 0) && (a > 0)) grad_event_terms[k] += 1 / a;
    k = evtid2base_probunits_param[evt_id];
    if ((k >= 0) && std::isfinite(b)) {
        double sigma = probunits2probability(b, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]),
               dsigma = probunits2probability_derivative(b, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
        if (sigma > 0) grad_event_terms[k] += dsigma / sigma;
    }
}

/** Add the rate of `count` pairs covered by the summary event of an event type to the total effective rate (and gradient).
 */
inline void add_summary_shares (int evt_id, double count)
{
    add_effective_rate(count * evtid2summary_single_er[evt_id]);
    register_summary_share_gradient(evt_id, count);
}

/** Subtract the rate of `count` pairs covered by the summary event of an event type from the total effective rate (and gradient).
 */
inline void subtract_summary_shares (int evt_id, double count)
{
    subtract_effective_rate(count * evtid2summary_single_er[evt_id]);
    register_summary_share_gradient(evt_id, -count);
}

/** Compute and register an event's effective rate and, in simulation mode, draw the tentative time
 *  at which it would happen if nothing changes in between and put it into the schedule.
 *
 *  In replay mode (scheduling_enabled == false), no random numbers are drawn and the schedule is not used;
 *  the event is only marked as registered by setting its time to infinity.
 */
inline void _schedule_event (
        event& ev,         ///< [in] the event to schedule
        event_data* evd_,  ///< [in] its data
        int evt_id         ///< [in] the dense id of its event type (see evt_id_of())
        )
{
    assert(evd_ == &ev2data.at(ev));
    assert(evt_id >= 0);
    rate ar = total_attempt_rate(evd_);
    if (ar < 0.0) throw "negative attempt rate";
    auto spu = total_success_probunits(evd_);
    double left_tail = evtid2left_tail[evt_id], right_tail = evtid2right_tail[evt_id], scale = evtid2scale[evt_id];
    bool is_summary = event_is_summary(ev);

    // compute the effective rate and register it in the total (and in the gradient of the total):
    rate er;
    if (is_summary)
    {
        // the summary event's effective rate uses the base success probability units
        // (the actual success of an attempt is tested later in pop_next_event):
        er = evd_->effective_rate = effective_rate(ar, spu, left_tail, right_tail, scale);
        assert (er < INFINITY);
        add_effective_rate(er);
        // NOTE: in update_adjacent_events, particular events must get update total effective rate properly as the difference to this!
        // (the gradient contributions of a summary event never change and are registered once in init_events)
    }
    else if (spu == -INFINITY)  // impossible event
    {
        er = evd_->effective_rate = 0;
    }
    else if (ar < INFINITY)
    {
        er = evd_->effective_rate = effective_rate(ar, spu, left_tail, right_tail, scale);
        assert (er < INFINITY);
        add_effective_rate(er);
        register_event_rate_gradient(evd_, evt_id, +1);
    }
    else  // event should happen "right away"
    {
        er = evd_->effective_rate = INFINITY;
        add_effective_rate(er);
    }

    if (!scheduling_enabled)
    {
        // replay mode: the event is registered but no tentative time is drawn:
        evd_->t = INFINITY;
        return;
    }

    // draw the tentative time at which the event would happen if nothing changes in between:
    timepoint t;
    if (is_summary)
    {
        // use a common upper bound to the actual effective rate for scheduling:
        t = current_t + exponential(random_variable) / (ar * evtid2summary_max_success_probability[evt_id]);
        if (verbose) cout << "         (re)scheduling " << ev << ": summary event, attempt rate " << ar << " → attempt at t=" << t << ", test success then" << endl;
    }
    else if (spu == -INFINITY)
    {
        t = INFINITY;
        if (debug) cout << "         (re)scheduling " << ev << ": zero success probability → t=" << t << endl;
    }
    else if (ar < INFINITY)
    {
        // draw time interval after which it would happen if nothing changes in between, and add it to the current time:
        t = current_t + exponential(random_variable) / er;
        if (verbose) {
            if (t == INFINITY) {
                if (debug) cout << "         (re)scheduling " << ev << ": zero effective rate → t=" << t << endl;
            }
            else cout << "         (re)scheduling " << ev << ": ar " << ar << ", spu " << spu << " → eff. rate " << er << " → next at t=" << t << endl;
        }
    }
    else
    {
        // to make sure that all "immediate" events occur in random order,
        // they are formally scheduled at some random "past" time:
        t = current_t - abs(1 + current_t) * uniform(random_variable);
        if (verbose) cout << "         (re)scheduling " << ev << ": ar inf, spu > 0 → eff. rate inf → next \"immediately\" at t=" << t << endl;
    }
    if (t == INFINITY)
    {
        // replace INFINITY by some unique finite but non-reached time point:
        t = never_t * (1 + uniform(random_variable));
    }
    // store time:
    evd_->t = t;
    t2ev[t] = ev;
}

inline void schedule_event (event& ev, event_data* evd_, int evt_id)
{
    assert(evd_ == &ev2data.at(ev));
    if (event_is_scheduled(ev, evd_)) throw "event already scheduled";
    assert(!event_is_scheduled(ev, evd_));
    _schedule_event(ev, evd_, evt_id);
    if (debug) verify_data_consistency();
}

/** Remove a scheduled event from the schedule and its rate (and gradient contributions) from the totals,
 *  but keep its data, so that the data can be modified and the event be scheduled anew (or erased).
 *
 *  Must be called BEFORE the event's attempt rate or success probunits are modified,
 *  since the gradient contributions must be removed with the same values they were registered with.
 */
inline void unschedule_event (event& ev, event_data* evd_, int evt_id)
{
    assert(evd_ == &ev2data.at(ev));
    assert(event_is_scheduled(ev, evd_));
    if (scheduling_enabled) t2ev.erase(evd_->t);
    bool is_summary = event_is_summary(ev);
    subtract_effective_rate(evd_->effective_rate, !is_summary);
    // (the gradient contributions of a summary event never change, so they are only registered once in init_events)
    if (!is_summary) register_event_rate_gradient(evd_, evt_id, -1);
    evd_->t = -INFINITY;
}

/** Reschedule an event whose data has NOT changed since it was scheduled (e.g. a summary event after an attempt).
 */
inline void reschedule_event (event& ev, event_data* evd_, int evt_id)
{
    unschedule_event(ev, evd_, evt_id);
    _schedule_event(ev, evd_, evt_id);
    if (debug) verify_data_consistency();
}

void add_event (event& ev);

void remove_event (event& ev, event_data* evd_);

void conditionally_remove_event(event& ev);

void update_adjacent_events (event& ev);

void add_reverse_event (event& old_ev);

void perform_event (event& ev, event_data* evd_);

void advance_time (timepoint t);

void finish_time ();

bool pop_next_event ();

#endif
