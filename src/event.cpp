/** Handling of events.
 *
 * An event is a tuple {ec, e1, rat13, e3} where ec is the event class
 * (EC_EST or EC_TERM for establishment or termination of a relationship,
 * EC_ACT for the occurrence of an act), e1 and e3 are the affected entities,
 * and rat13 is the relationship or action type.
 *
 * An event type is a tuple (ec, et1, rat13, et3} where et1, et3 are two entity types.
 *
 * An event's variable data (current attempt rate, success probability units, and next scheduled timepoint)
 * is stored in a separate struct of type event_data.
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Mar 30, 2020
 *
 * @file
 */

#include <assert.h>
#include <iostream>

#include "global_variables.h"
#include "entity.h"
#include "angle.h"
#include "link.h"
#include "probability.h"
#include "event.h"
#include "io.h"
#include "debugging.h"

/** Return whether a future instance of the event is scheduled.
 */
bool event_is_scheduled (
        event& ev, ///< [in] the event to check, passed by reference for performance
        event_data* evd_ ///< [in] the corresponding variable event data, passed by pointer for performance
        )
{
    assert(evd_ == &ev2data.at(ev));
    return (evd_->t > -INFINITY);
}

void _schedule_event (event& ev, event_data* evd_, double left_tail, double right_tail)
{
    assert(evd_ == &ev2data.at(ev));
    timepoint t;
    auto ar = evd_->attempt_rate;
    if (ar < 0.0) throw "negative attempt rate";
    assert (ar >= 0.0);
    if (verbose) cout << "         (re)scheduling " << ev << ": ";
    if (event_is_summary(ev)) { // summary event: use only attempt rate (success will be tested in pop_next_event):
        t = current_t + exponential(random_variable) / (ar * ev2max_success_probability[ev]);
        if (verbose) cout << "summary event, attempt rate " << ar << " → attempt at t=" << t << ", test success then" << endl;
    } else { // particular event: use effective rate:
        auto spu = evd_->success_probunits;
        if (spu == -INFINITY) {
            t = INFINITY;
            if (verbose) cout << "zero success probability → t=" << t << endl;
        } else if (ar < INFINITY) {
            auto er = effective_rate(ar, spu, left_tail, right_tail);
            t = current_t + exponential(random_variable) / er;
            if (verbose) cout << "ar " << ar << ", spu " << spu << " → eff. rate " << er << " → next at t=" << t << ((t==INFINITY) ? " (never)" : "") << endl;
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

void schedule_event (event& ev, event_data* evd_, double left_tail, double right_tail)
{
    assert(evd_ == &ev2data.at(ev));
    if (event_is_scheduled(ev, evd_)) throw "event already scheduled";
    assert(!event_is_scheduled(ev, evd_));
    _schedule_event(ev, evd_, left_tail, right_tail);
    if (debug) verify_data_consistency();
}

void reschedule_event (event& ev, event_data* evd_, double left_tail, double right_tail)
{
    assert(evd_ == &ev2data.at(ev));
    assert(event_is_scheduled(ev, evd_));
    t2ev.erase(evd_->t);
    _schedule_event(ev, evd_, left_tail, right_tail);
    if (debug) verify_data_consistency();
}

void add_event (event& ev)
{
    auto ec = ev.ec; auto e1 = ev.e1, e3 = ev.e3; auto rat13 = ev.rat13;
    assert ((rat13 != RT_ID) && (e1 != e3));
    auto et1 = e2et[e1], et3 = e2et[e3];
    event_type evt = { .ec = ec, .et1 = et1, .rat13 = rat13, .et3 = et3 };

    if (possible_evts.count(evt) > 0) { // event can happen at all:

        if (debug) cout << "     adding event: " << ev << endl;
        // find and store attempt rate and success probunit by looping through adjacent legs and angles:
        rate ar = evt2base_attempt_rate[evt];
        auto spu = evt2base_probunits[evt];
        auto outs1 = e2outs[e1];
        auto ins3 = e2ins[e3];
        // legs:
        if (ec == EC_TERM) {
            for (auto& [rat12, e2] : outs1) {
                influence_type inflt = { .evt = evt, .at = { .rat12 = rat12, .et2 = e2et[e2], .rat23 = NO_RAT } };
                ar += _inflt2attempt_rate[INFLT(inflt)];
                spu += _inflt2delta_probunits[INFLT(inflt)];
            }
            for (auto& [e2, rat23] : ins3) {
                influence_type inflt = { .evt = evt, .at = { .rat12 = NO_RAT, .et2 = e2et[e2], .rat23 = rat23 } };
                ar += _inflt2attempt_rate[INFLT(inflt)];
                spu += _inflt2delta_probunits[INFLT(inflt)];
            }
        }

        if (debug) {
            cout << "outs:" << endl;
            for (auto& [rat13, e3] : outs1) cout << " " << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << endl;
            cout << "ins:" << endl;
            for (auto& [e1, rat13] : ins3) cout << " " << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << endl;
        }

        // angles:
        int na = 0; // number of influencing angles
        angles as = leg_intersection(e1, outs1, ins3, e3);
        for (auto a_it = as.begin(); a_it != as.end(); a_it++) {
            influence_type inflt = {
                    .evt = evt,
                    .at = { .rat12 = a_it->rat12, .et2 = e2et[a_it->e2], .rat23 = a_it->rat23 }
            };
            if (debug) cout << "      influences of angle \"" << e2label[e1] << " " << rat2label[a_it->rat12] << " " << e2label[a_it->e2] << " " << rat2label[a_it->rat23] << " " << e2label[e3] << "\":" << endl;
            auto dar = _inflt2attempt_rate[INFLT(inflt)];
            auto dspu = _inflt2delta_probunits[INFLT(inflt)];
            if (COUNT_ALL_ANGLES || (dar != 0.0) || (dspu != 0.0)) { // angle can influence event
                na++;
                if (debug) {
                    if (dar != 0.0) cout << "       on attempt rate:" << dar << endl;
                    if (dspu != 0.0) cout << "       on success probunit:" << dspu << endl;
                }
                ar += dar;
                spu += dspu;
            } else if (debug) cout << "       none" << endl;
        }

        // only add a non-termination event individually if at least one angle existed
        // (pure spontaneous non-termination events are handled summarily via summary events to keep maps sparse):
        if ((ec == EC_TERM) || (na > 0)) {
            assert (ev2data.count(ev) == 0);
            ev2data[ev] = { .n_angles = na, .attempt_rate = max(0.0, ar), .success_probunits = spu, .t = -INFINITY };
            if (debug) cout << "      attempt rate " << ar << ", success prob. " << probunits2probability(spu, evt2left_tail.at(evt), evt2right_tail.at(evt)) << endl;
            schedule_event(ev, &ev2data[ev], evt2left_tail.at(evt), evt2right_tail.at(evt));
            if (debug) { verify_data_consistency(); verify_angle_consistency(); }
        } else {
            if (debug) cout << "      covered by summary event, not scheduled separately" << endl;
        }

    } else {
        if (debug) cout << "     not adding impossible event: " << ev << endl;
    }
}

void remove_event (event& ev, event_data* evd_)
{
    assert (event_is_scheduled(ev, evd_));
    t2ev.erase(evd_->t);
    ev2data.erase(ev);
    if (debug) cout << "        removed event: " << ev << " scheduled at " << evd_->t << endl;
}

void conditionally_remove_event(event& ev)
{
    if (ev2data.count(ev) == 1) {
        auto evd_ = &ev2data.at(ev);
        remove_event(ev, evd_);
    }
}

#define IF_12 if (which == 0)
#define IF_23 if (which == 1)

void update_adjacent_events (event& ev)
{
    auto ec_ab = ev.ec; auto ea = ev.e1, eb = ev.e3; auto rab = ev.rat13;
    entity e1, e2, e3; entity_type et1, et2, et3; relationship_or_action_type rat12, rat23;
    if (debug) cout << "  updating adjacent events of " << ev << endl;
    // loop through all adjacent events e1->e3:

    if (debug) cout << "   angles with this as 1st leg:" << endl;
    e1 = ea; rat12 = rab; e2 = eb;
    et1 = e2et[e1]; et2 = e2et[e2];
    auto outlegs = e2outs[eb]; // rat23, e3
    for (auto& l : outlegs) {
        rat23 = l.rat_out; e3 = l.e_target; et3 = e2et[e3];
        if (e1 != e3) { // since we allow no self-links except identity
            add_or_delete_angle(ec_ab, e1, et1, rat12, e2, et2, rat23, e3, et3);
        }
    }

    if (debug) cout << "   angles with this as 2nd leg:" << endl;
    auto inlegs = e2ins[ea]; // e1, rat12
    e2 = ea; rat23 = rab; e3 = eb;
    et2 = e2et[e2]; et3 = e2et[e3];
    for (auto& l : inlegs) {
        e1 = l.e_source; rat12 = l.rat_in; et1 = e2et[e1];
        if (e1 != e3) { // since we allow no self-links except identity
            add_or_delete_angle(ec_ab, e1, et1, rat12, e2, et2, rat23, e3, et3);
        }
    }
}

void add_reverse_event (event& old_ev)
{
    auto e1 = old_ev.e1, e3 = old_ev.e3; auto rat13 = old_ev.rat13;
    event_class ec = (old_ev.ec == EC_TERM) ? EC_EST : EC_TERM;
    event ev = { .ec = ec, .e1 = e1, .rat13 = rat13, .e3 = e3 };
    if (debug) cout << "    adding reverse event: " << ev << endl;
    add_event(ev);
}

void perform_event (event& ev)
{
    if (debug) cout << " performing event: " << ev << endl;

    auto ec = ev.ec; auto e1 = ev.e1, e3 = ev.e3; auto rat13 = ev.rat13, rat31 = rat2inv.at(rat13);
    assert (rat13 != RT_ID);
    link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 };

    // FIRST add the reverse event, so that its n_angles will reflect the situation before the change:
    add_reverse_event(ev);
    // THEN add or remove the link to perform the change:
    if (ec == EC_EST) {
        add_link(l);
    } else {
        del_link(l);
    }
    // FINALLY update all adjacent events (including the reverse event) to reflect the change:
    update_adjacent_events(ev);

    // also perform companion event that affects inverse link:
    if (rat31 != NO_RAT) {
        event companion_ev = { .ec = ec, .e1 = e3, .rat13 = rat31, .e3 = e1 };
        link inv_l = { .e1 = e3, .rat13 = rat31, .e3 = e1 }; // inverse link
        if (ev2data.count(companion_ev) == 1) {
            if (debug) cout << " unscheduling companion event: " << companion_ev << endl;
            auto companion_evd_ = &ev2data.at(companion_ev);
            remove_event(companion_ev, companion_evd_);
        }
        // FIRST add the reverse event, so that its n_angles will reflect the situation before the change:
        add_reverse_event(companion_ev);
        // THEN add or remove the link to perform the change:
        if (ec == EC_EST) {
            if (debug) cout << " performing companion event: adding inverse link \"" << e2label[e3]
                << " " << rat2label[rat31] << " " << e2label[e1] << "\"" << endl;
            add_link(inv_l);
        } else {
            if (debug) cout << " performing companion event: deleting inverse link \"" << e2label[e3]
                << " " << rat2label[rat31] << " " << e2label[e1] << "\"" << endl;
            del_link(inv_l);
        }
        // FINALLY update all adjacent events (including the reverse event) to reflect the change:
        update_adjacent_events(companion_ev);
    }
    if (debug) {
        dump_links();
        verify_angle_consistency();
    }
}

bool pop_next_event ()
{
    bool found = false;
    // find next event:
    while ((!found) && (current_t < max_t)) {

        // get handle of earliest next scheduled event:
        auto tev = t2ev.begin();
        if (tev == t2ev.end()) { // no events are scheduled --> model has converged
            log_status();
            // jump to end
            current_t = max_t;
            return false;
        }

        // get corresponding timepoint:
        timepoint t = tev->first;
        if (t >= max_t) { // no events before max_t are scheduled
            if (!quiet) {
                if (t < INFINITY) cout << "next event would happen after time limit at t=" << t << endl;
                else cout << "no further events are scheduled." << endl;
            }
            // jump to end:
            current_t = max_t;
            return false;
        }

        // get the corresponding event:
        event ev = tev->second;
        if (t > current_t) {
            // advance model time to time of event:
            current_t = t;
        } // otherwise it's an event happening "right now" scheduled formally for a past time to ensure a random order of those events.

        if (event_is_summary(ev)) { // event is the attempt of a summary event, so has only types specified

            event summary_ev = ev;
            entity_type et1 = summary_et1(summary_ev), et3 = summary_et3(summary_ev);
            if (debug) cout << "at t=" << current_t << " summary event " << summary_ev << " :" << endl;

            // draw actual entities at random from given types:
            auto e1 = random_entity(et1), e3 = random_entity(et3);
            auto rat13 = summary_ev.rat13;
            event_type evt = { .ec = EC_EST, .et1 = et1, .rat13 = rat13, .et3 = et3 };
            link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 };

            if (link_exists(l)) {
                if (verbose) cout << "at t=" << current_t << ", link to establish \"" << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << "\" existed already" << endl;
            } else if (e1 == e3) {
                if (verbose) cout << "at t=" << current_t << ", entities to link were equal and are thus not linked" << endl;
            } else { // link can be established
                event actual_ev = { .ec = EC_EST, .e1 = e1, .rat13 = rat13, .e3 = e3 };
                if (ev2data.count(actual_ev) > 0) {
                    // the event was scheduled separately, e.g. due to an angle, so we do not perform it now.
                    if (verbose) cout << "at t=" << current_t << " " << actual_ev << " is scheduled separately at t=" << ev2data.at(actual_ev).t << ", so not performed now." << endl;
                } else { // check if attempt succeeds
                    // compile success units:
                    auto spu = evt2base_probunits.at(evt);
                    for (auto& [rat12, e2] : e2outs[e1]) {
                        influence_type inflt = { .evt = evt, .at = { .rat12 = rat12, .et2 = e2et[e2], .rat23 = NO_RAT } };
                        if (inflt2delta_probunits.count(inflt) > 0) spu += inflt2delta_probunits.at(inflt);
                    }
                    for (auto& [e2, rat23] : e2ins[e3]) {
                        influence_type inflt = { .evt = evt, .at = { .rat12 = NO_RAT, .et2 = e2et[e2], .rat23 = rat23 } };
                        if (inflt2delta_probunits.count(inflt) > 0) spu += inflt2delta_probunits.at(inflt);
                    }
                    // since the scheduling rate already contained the factor ev2max_sp[ev],
                    // we need to divide the success probability by it here:
                    probability conditional_success_probability =
                            probunits2probability(spu, evt2left_tail.at(evt), evt2right_tail.at(evt)) / ev2max_success_probability[ev];
                    // check if success:
                    if (uniform(random_variable) < conditional_success_probability) { // success
                        // "return" event:
                        current_ev = actual_ev;
                        if (!quiet) log_status();
                        found = true;
                    } else {
                        if (verbose) cout << "at t=" << current_t << " " << actual_ev << " did not succeed" << endl;
                    }
                }
            }
            // set next_occurrence of this summary event:
            reschedule_event(summary_ev, &ev2data.at(summary_ev), evt2left_tail.at(evt), evt2right_tail.at(evt));

        } else { // event has specific entities

            // "return" event:
            current_ev = ev;
            if (!quiet) log_status();
            auto evd_ = current_evd_ = &ev2data.at(ev);
            // remove it from all relevant data:
            remove_event(ev, evd_);
            found = true;
        }
    }
    return found;
}



