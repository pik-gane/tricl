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
#include "legs.h"
#include "actions.h"

/** Perform all necessary changes in state and event data
 *  due to the addition or deletion of an angle, or a change of its weight.
 *
 *  Iterate through all events that might be influenced by this angle,
 *  update their attempt rates and success probability units (or, for angles involving actions,
 *  their decaying influence weights, see actions.h), and (re)schedule them based on the new rates.
 *
 *  This is one of the performance bottleneck functions
 *  since it is called many times by \ref update_adjacent_angles().
 */
inline void add_or_delete_angle (
        int dcount,                         ///< [in] +1 if the angle is added, -1 if it is deleted, 0 if only its weight changes
        double dw,                          ///< [in] change of the weight of the angle (only used for angles involving actions)
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
    if (debug) cout << "    angle change " << dcount << " (weight " << dw << ") \"" << e2label[e1] << " " << rat2label[rat12] << " "
            << e2label[e2] << " " << rat2label[rat23] << " " << e2label[e3] << "\"" << endl;

    // update total no. of angles (if non-id.):
    if (!((e1 == e2) || (e2 == e3) || (e3 == e1))) n_angles += dcount;

    // iterate through all possible source-target relationship or action types:
    for (auto& rat13 : ets2rats[(int) et1 * n_et_slots + (int) et3])
    {
        // construct the type of the corresponding event whose data might need an update:
        event_class ec13;
        if (any_actions && rat_is_action[rat13]) ec13 = EC_ACT;
        else ec13 = (e2outs[e1].count({ .rat_out = rat13, .e_target = e3 }) > 0) ? EC_TERM : EC_EST;
        if (debug) cout << "     possibly updating event: " << ec2label[ec13] <<  " \"" << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << "\"" << endl;

        // only continue if the event type can happen at all:
        int evt_id = evt_id_of(ec13, et1, rat13, et3);
        if (evt_id < 0) continue;
        event ev = { .ec=ec13, e1, rat13, e3 };

        // get influence of angle on event, and only continue if influence is nonzero:
        int idx = inflt_index(evt_id, rat12, et2, rat23);
        int j = evtid_at2infl[idx];
        if (j < 0) {
            if (debug) cout << "       but angle may not influence attempt or success" << endl;
            continue;
        }
        bool decaying = any_actions && (evtid2infl_decay[evt_id][j] > 0);
        if (!decaying && (dcount == 0)) continue;
        auto dar = inflt_attempt_rate[idx];
        auto dspu = inflt_delta_probunits[idx];
        if (debug) cout << "       angle may influence attempt or success" << endl;

        auto it = ev2data.find(ev);
        if (it == ev2data.end())  // event is not scheduled individually
        {
            if (dcount <= 0) continue;  // (an angle of a pair covered by a summary event cannot be removed or reweighted)
            if (debug) cout << "        event will be scheduled newly" << endl;
            if (ec13 != EC_TERM)  // event is ALSO covered by a summary event
            {
                // subtract that part covered by the summary event from the total effective rate:
                if (evtid2summary[evt_id]) subtract_pair_share(evt_id, e1, e3);
            }
            auto evd_ = &ev2data[ev];  // generates a new (zero-initialised) event_data object
            index_event(ev);
            evd_->n_angles = 1;
            evd_->n_infl[j] = 1;
            add_attempt_contribution(evd_, evtid2base_attempt_rate[evt_id]);
            add_probunits_contribution(evd_, evtid2base_probunits[evt_id]);
            if (decaying) {
                auto& tv = ev2tv[ev];
                tv = {};
                tv.t_ref = current_t;
                tv.W[j] = max(0.0, dw);
            } else {
                add_attempt_contribution(evd_, dar);
                add_probunits_contribution(evd_, dspu);
            }
            if (any_leg_influences) add_leg_contributions(evt_id, e1, rat13, e3, *evd_);
            schedule_event(ev, evd_, evt_id);
            continue;
        }

        if (debug) cout << "        event will be rescheduled" << endl;
        auto evd_ = &(it->second);
        // unregister the event's rate and gradient contributions BEFORE modifying its data
        // (for time-varying events, this also accumulates the integral of the rate up to now):
        unschedule_event(ev, evd_, evt_id);
        if (dcount < 0) {
            assert (evd_->n_angles > 0);  // since angle must have been added earlier to be removed now
            assert (evd_->n_infl[j] > 0);
        }
        evd_->n_angles += dcount;
        evd_->n_infl[j] += dcount;
        if (decaying) {
            auto tit = ev2tv.find(ev);
            if (tit == ev2tv.end()) {
                tit = ev2tv.emplace(ev, tv_data{}).first;
                tit->second.t_ref = current_t;
            }
            auto& tv = tit->second;
            tv.W[j] = (evd_->n_infl[j] == 0) ? 0.0 : max(0.0, tv.W[j] + dw);
            // drop the time-varying data if no decaying influence is left:
            bool any_left = false;
            const auto& dec = evtid2infl_decay[evt_id];
            for (size_t k = 0; k < dec.size(); k++) if ((dec[k] > 0) && (evd_->n_infl[k] > 0)) any_left = true;
            if (!any_left) ev2tv.erase(tit);
        } else if (dcount > 0) {
            add_attempt_contribution(evd_, dar);
            add_probunits_contribution(evd_, dspu);
        } else {
            remove_attempt_contribution(evd_, dar);
            remove_probunits_contribution(evd_, dspu);
        }
        if ((ec13 != EC_TERM) && (evd_->n_angles == 0))  // only spontaneous non-termination event is left:
        {
            if (debug) cout << "        last angle was removed, so event will be removed because it is covered by a summary event" << endl;
            // remove specific event:
            unindex_event(ev);
            if (any_actions) ev2tv.erase(ev);
            ev2data.erase(it);
            // add that part covered by the summary event to the total effective rate:
            if (evtid2summary[evt_id]) add_pair_share(evt_id, e1, e3);
        }
        else
        {
            _schedule_event(ev, evd_, evt_id);
            if (debug) verify_data_consistency();
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
