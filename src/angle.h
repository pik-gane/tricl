// make sure this file is only included once:
#ifndef INC_ANGLE_H
#define INC_ANGLE_H

/** Performance-critical inline functions for handling of angles.
 *
 *  \file
 *
 *  See \ref data_model.h for how a \ref tricl::angle relates to other tricl datatypes.
 */

#include "global_variables.h"
#include "probability.h"
#include "debugging.h"
#include "event.h"
#include "io.h"

/** Perform all necessary changes in state and event data
 *  due to the addition or deletion of an angle.
 *
 *  Iterate through all events that might be influenced by this angle,
 *  update their attempt rates and success probability units,
 *  and (re)schedule them based on the new rates.
 *
 *  This is one of the performance bottleneck functions
 *  since it is called many times by \ref update_adjacent_events().
 */
inline void add_or_delete_angle (
        event_class ec_angle,               ///< [in] class of the event that caused the change: EC_EST adds an angle, EC_TERM deletes one
        entity e1,                          ///< [in] source entity
        entity_type et1,                    ///< [in] its type
        relationship_or_action_type rat12,  ///< [in] source-to-middle relationship or action type
        entity e2,                          ///< [in] middle entity
        entity_type et2,                    ///< [in] its type
        relationship_or_action_type rat23,  ///< [in] middle-to-target relationship or action type
        entity e3,                          ///< [in] target entity
        entity_type et3                     ///< [in] its type
        )
{
    if (debug) cout << "    " << ec2label[ec_angle] << " \"" << e2label[e1] << " " << rat2label[rat12] << " "
            << e2label[e2] << " " << rat2label[rat23] << " " << e2label[e3] << "\"" << endl;

    // update total no. of angles (if non-id.):
    n_angles += ((e1 == e2) || (e2 == e3) || (e3 == e1)) ? 0 : (ec_angle == EC_EST) ? 1 : -1;

    // iterate through all possible source-target relationship or action types:
    for (auto& rat13 : ets2rats[(int) et1 * n_et_slots + (int) et3])
    {
        bool link13_exists = (e2outs[e1].count({ .rat_out = rat13, .e_target = e3 }) > 0);

        // construct the type of the corresponding event whose data might need an update:
        event_class ec13 = link13_exists ? EC_TERM : EC_EST;
        if (debug) cout << "     possibly updating event: " << ec2label[ec13] <<  " \"" << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << "\"" << endl;

        // only continue if the event type can happen at all:
        int evt_id = evt_id_of(ec13, et1, rat13, et3);
        if (evt_id >= 0)
        {
            event ev = { .ec=ec13, e1, rat13, e3 };
            if (debug) cout << "      event type " << evtid2evt[evt_id] << " has a base success prob. of "
                    << probunits2probability(evtid2base_probunits[evt_id], evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]) << endl;

            // get influence of angle on event:
            int idx = inflt_index(evt_id, rat12, et2, rat23);
            auto dar = inflt_attempt_rate[idx];
            auto dspu = inflt_delta_probunits[idx];

            // only continue if influence is nonzero:
            if (COUNT_ALL_ANGLES || (dar != 0.0) || (dspu != 0.0))
            {
                if (debug) cout << "       angle may influence attempt or success" << endl;
                if (ec_angle == EC_EST)  // angle is added:
                {
                    if (ev2data.count(ev) == 0)  // event is not already scheduled
                    {
                        if (debug) cout << "        event will be scheduled newly" << endl;
                        if (ec13 != EC_TERM)  // event is ALSO covered by a summary event
                        {
                            // subtract that part covered by the summary event from the total effective rate:
                            subtract_effective_rate(evtid2summary_single_er[evt_id]);
                        }
                        auto evd_ = &ev2data[ev];  // generates a new (zero-initialised) event_data object
                        evd_->n_angles = 1;
                        add_attempt_contribution(evd_, evtid2base_attempt_rate[evt_id]);
                        add_attempt_contribution(evd_, dar);
                        add_probunits_contribution(evd_, evtid2base_probunits[evt_id]);
                        add_probunits_contribution(evd_, dspu);
                        schedule_event(ev, evd_, evt_id);
                    }
                    else  // event is already scheduled
                    {
                        if (debug) cout << "        event will be rescheduled" << endl;
                        auto evd_ = &(ev2data.at(ev));
                        evd_->n_angles += 1;
                        add_attempt_contribution(evd_, dar);
                        add_probunits_contribution(evd_, dspu);
                        reschedule_event(ev, evd_, evt_id);
                    }
                }
                else  // angle is removed
                {
                    auto evd_ = &(ev2data.at(ev));
                    assert (evd_->n_angles > 0);  // since angle must have been added earlier to be removed now
                    evd_->n_angles -= 1;
                    remove_attempt_contribution(evd_, dar);
                    remove_probunits_contribution(evd_, dspu);
                    if ((ec13 != EC_TERM) && (evd_->n_angles == 0))  // only spontaneous non-termination event is left:
                    {
                        if (debug) cout << "        last angle was removed, so event will be removed because it is covered by a summary event" << endl;
                        // remove specific event:
                        remove_event(ev, evd_);  // event must have been scheduled earlier when angle was added
                        // add that part covered by the summary event to the total effective rate:
                        add_effective_rate(evtid2summary_single_er[evt_id]);
                    }
                    else
                    {
                        if (debug) cout << "        event will be rescheduled" << endl;
                        reschedule_event(ev, evd_, evt_id);  // event must have been scheduled earlier when angle was added
                    }
                }
            }
            else
            {
                if (debug) cout << "       but angle may not influence attempt or success" << endl;
            }
        }
    }
}

/** Advance an iterator over a range of legs sorted by their "other" entity
 *  to the first leg whose other entity is >= key, using galloping (exponential + binary) search.
 *
 *  Precondition: it != end and e2of(*it) < key.
 *
 *  This makes the intersection of a small leg set with a large one (e.g. of a hub entity
 *  that encodes a state and is linked to all agents) cost O(small * log(large)) instead of O(large).
 */
template <class It, class E2Of>
inline It advance_legs_to (It it, It end, entity key, E2Of e2of)
{
    size_t n = end - it, step = 1, pos = 0;  // invariant: e2of(it[pos]) < key
    while ((pos + step < n) && (e2of(it[pos + step]) < key)) {
        pos += step;
        step *= 2;
    }
    // now the answer lies in (pos, min(pos + step, n)]:
    size_t lo = pos + 1, hi = min(pos + step, n);
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (e2of(it[mid]) < key) lo = mid + 1; else hi = mid;
    }
    return it + lo;
}

/** Compare each \ref outleg of e1 with each \ref inleg of e3 to find each \ref angle from e1 to e3.
 *
 *  This is one of the performance bottleneck functions
 *  since it is called by \ref add_event() for every event that is added.
 *
 *  Both leg sets are sorted by the middle entity e2 (see the comparison operators of \ref inleg and \ref outleg),
 *  so the angles can be found by a merge-like intersection, in which the iterator over the larger set
 *  is advanced by galloping search (see \ref advance_legs_to()).
 *
 *  The found angles are appended to the caller-provided buffer (which is cleared first) in the order of
 *  increasing e2, and for each e2 ordered by inleg first and outleg second.
 *  (This order determines the summation order of floating point contributions and must not be changed
 *  without regenerating the regression test references.)
 */
inline void get_angles (
        const entity e1,         ///< [in] source entity
        const outleg_set& out1,  ///< [in] set of outlegs of source entity
        const inleg_set& in3,    ///< [in] set of inlegs of target entity
        const entity e3,         ///< [in] target entity
        angle_vec& result        ///< [out] buffer to store the found angles in (reused between calls to avoid allocations)
        )
{
    result.clear();
    auto o = out1.begin(), oend = out1.end();
    auto i = in3.begin(), iend = in3.end();
    auto oe2 = [](const outleg& l) { return l.e_target; };
    auto ie2 = [](const inleg& l) { return l.e_source; };
    while ((o != oend) && (i != iend))
    {
        entity eo = o->e_target, ei = i->e_source;
        if (eo < ei)
        {
            o = advance_legs_to(o, oend, ei, oe2);
        }
        else if (ei < eo)
        {
            i = advance_legs_to(i, iend, eo, ie2);
        }
        else  // common middle entity e2 = eo = ei
        {
            // find the ends of the blocks of legs with this middle entity (at most n_rats legs each):
            auto oblockend = o;
            while ((oblockend != oend) && (oblockend->e_target == eo)) ++oblockend;
            auto iblockend = i;
            while ((iblockend != iend) && (iblockend->e_source == eo)) ++iblockend;
            // store all combinations:
            for (auto ii = i; ii != iblockend; ++ii) {
                for (auto oo = o; oo != oblockend; ++oo) {
                    result.push_back({ .rat12 = oo->rat_out, .e2 = eo, .rat23 = ii->rat_in });
                    if (debug) cout << "      angle: " << e2label[e1] << " " << rat2label[oo->rat_out] << " "
                          << e2label[eo] << " " << rat2label[ii->rat_in] << " " << e2label[e3] << endl;
                }
            }
            o = oblockend;
            i = iblockend;
        }
    }
}

#endif
