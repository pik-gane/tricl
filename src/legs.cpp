/** Influences of legs on events (see legs.h).
 *
 *  \file
 */

#include "global_variables.h"
#include "probability.h"
#include "entity.h"
#include "event.h"
#include "legs.h"

bool track_legs = false;
vector<legkey_info> legkeys;
vector<int> slot2legkey;
vector<int> e2legdeg;
vector<vector<evt_leg>> evtid2outlegs, evtid2inlegs;
vector<char> evtid_is_leg_summary;
vector<leg_summary_data> evtid2lsd;
unordered_map<link_type, link_sampler> lt2sampler;
vector<unordered_set<event>> e2src_evs, e2tgt_evs;
vector<int> evtid2summary;
static vector<event> evtid2summary_ev;           // the summary event of each leg summary event type
static vector<vector<int>> legkey2summaries;     // leg summary event types influenced by each leg type
static vector<char> evtid_dirty;                 // leg summary event types whose summary event must be rescheduled
static vector<int> dirty_list;

/** Add the gradient of the class rate lambda(u, v) times a weight to \ref grad_rate. */
void register_class_gradient (int evt_id, const vector<int>& u, const vector<int>& v, double weight)
{
    if (!compute_gradient || (weight == 0.0)) return;
    rate A; probunits B;
    class_rate(evt_id, u, v, A, B);
    if (B == -INFINITY) return;
    double sigma = probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]),
           dsigma = probunits2probability_derivative(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
    int k = evtid2base_attempt_param[evt_id];
    if (k >= 0) grad_rate[k] += weight * sigma;
    k = evtid2base_probunits_param[evt_id];
    if ((k >= 0) && (dsigma != 0.0)) grad_rate[k] += weight * A * dsigma;
    const auto& ia = evtid2infl_attempt_param[evt_id];
    const auto& ib = evtid2infl_probunits_param[evt_id];
    auto legs = [&](const vector<evt_leg>& ls, const vector<int>& n) {
        for (size_t i = 0; i < ls.size(); i++) if (n[i] > 0) {
            int j = ls[i].j;
            if (ia[j] >= 0) grad_rate[ia[j]] += weight * n[i] * sigma;
            if ((ib[j] >= 0) && (dsigma != 0.0)) grad_rate[ib[j]] += weight * n[i] * A * dsigma;
        }
    };
    legs(evtid2outlegs[evt_id], u);
    legs(evtid2inlegs[evt_id], v);
}

/** Add the gradient of log lambda(u, v) to \ref grad_event_terms (for an event of a covered pair that happens now). */
void add_class_log_gradient (int evt_id, const vector<int>& u, const vector<int>& v)
{
    if (!compute_gradient) return;
    rate A; probunits B;
    class_rate(evt_id, u, v, A, B);
    if (!(A > 0) || (B == -INFINITY)) return;
    double sigma = probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]),
           dsigma = probunits2probability_derivative(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
    double dlogsigma = (sigma > 0) ? dsigma / sigma : 0.0;
    int k = evtid2base_attempt_param[evt_id];
    if (k >= 0) grad_event_terms[k] += 1 / A;
    k = evtid2base_probunits_param[evt_id];
    if (k >= 0) grad_event_terms[k] += dlogsigma;
    const auto& ia = evtid2infl_attempt_param[evt_id];
    const auto& ib = evtid2infl_probunits_param[evt_id];
    auto legs = [&](const vector<evt_leg>& ls, const vector<int>& n) {
        for (size_t i = 0; i < ls.size(); i++) if (n[i] > 0) {
            int j = ls[i].j;
            if (ia[j] >= 0) grad_event_terms[ia[j]] += n[i] / A;
            if (ib[j] >= 0) grad_event_terms[ib[j]] += n[i] * dlogsigma;
        }
    };
    legs(evtid2outlegs[evt_id], u);
    legs(evtid2inlegs[evt_id], v);
}

/** Add the rate (and gradient) of a pair that becomes covered by the summary event of its event type to the totals. */
void add_pair_share (int evt_id, entity e1, entity e3)
{
    if (!evtid_is_leg_summary[evt_id]) {
        add_summary_shares(evt_id, 1);
        return;
    }
    auto u = outdegs(evt_id, e1), v = indegs(evt_id, e3);
    rate A; probunits B;
    total_finite_effective_rate += class_rate(evt_id, u, v, A, B);
    register_class_gradient(evt_id, u, v, +1);
}

/** Subtract the rate (and gradient) of a pair that ceases to be covered by the summary event of its event type from the totals. */
void subtract_pair_share (int evt_id, entity e1, entity e3)
{
    if (!evtid_is_leg_summary[evt_id]) {
        subtract_summary_shares(evt_id, 1);
        return;
    }
    auto u = outdegs(evt_id, e1), v = indegs(evt_id, e3);
    rate A; probunits B;
    total_finite_effective_rate -= class_rate(evt_id, u, v, A, B);
    register_class_gradient(evt_id, u, v, -1);
}

/** Add the gradient of the log rate of an event of a covered pair that happens now. */
void add_covered_event_log_gradient (int evt_id, entity e1, entity e3)
{
    if (!evtid_is_leg_summary[evt_id]) add_summary_event_log_gradient(evt_id);
    else add_class_log_gradient(evt_id, outdegs(evt_id, e1), indegs(evt_id, e3));
}

/** Register an individually scheduled event in the per-entity indices. */
void index_event (const event& ev)
{
    if (!track_legs || event_is_summary(ev)) return;
    e2src_evs[ev.e1].insert(ev);
    e2tgt_evs[ev.e3].insert(ev);
}

/** Remove an individually scheduled event from the per-entity indices. */
void unindex_event (const event& ev)
{
    if (!track_legs || event_is_summary(ev)) return;
    e2src_evs[ev.e1].erase(ev);
    e2tgt_evs[ev.e3].erase(ev);
}

/** \returns the link type sampled by a leg channel of an event type. */
static link_type channel_link_type (const event_type& evt, const legkey_info& lk)
{
    return lk.out ? link_type{ evt.et1, lk.rat, lk.et2 } : link_type{ lk.et2, lk.rat, evt.et3 };
}

/** Set up the leg data (after init_types(), before init_events()). */
void init_legs ()
{
    legkeys.clear();
    slot2legkey.assign(n_at_slots, -1);
    evtid2outlegs.assign(n_evt_ids, {});
    evtid2inlegs.assign(n_evt_ids, {});
    for (int id = 0; id < n_evt_ids; id++) {
        for (size_t j = 0; j < evtid2infl_slots[id].size(); j++) {
            int slot = evtid2infl_slots[id][j];
            auto rat12 = (relationship_or_action_type) (slot / (n_et_slots * n_rat_slots));
            auto et2 = (entity_type) ((slot / n_rat_slots) % n_et_slots);
            auto rat23 = (relationship_or_action_type) (slot % n_rat_slots);
            if ((rat12 != NO_RAT) && (rat23 != NO_RAT)) continue;  // an angle
            if (slot2legkey[slot] < 0) {
                slot2legkey[slot] = legkeys.size();
                legkeys.push_back({ rat23 == NO_RAT, (rat23 == NO_RAT) ? rat12 : rat23, et2 });
            }
            evt_leg el = { (int) j, slot2legkey[slot], slot };
            if (rat23 == NO_RAT) evtid2outlegs[id].push_back(el); else evtid2inlegs[id].push_back(el);
        }
    }
    track_legs = !legkeys.empty();
    e2legdeg.assign(track_legs ? ((size_t) max_e + 1) * legkeys.size() : 0, 0);
    e2src_evs.assign(track_legs ? max_e + 1 : 0, {});
    e2tgt_evs.assign(track_legs ? max_e + 1 : 0, {});

    // which event types have summary events (establishment and act events with a positive basic or leg attempt rate):
    evtid2summary.assign(n_evt_ids, 0);
    evtid_is_leg_summary.assign(n_evt_ids, 0);
    evtid2lsd.assign(n_evt_ids, {});
    evtid2summary_ev.assign(n_evt_ids, {});
    evtid_dirty.assign(n_evt_ids, 0);
    legkey2summaries.assign(legkeys.size(), {});
    lt2sampler.clear();
    for (int id = 0; id < n_evt_ids; id++) {
        auto& evt = evtid2evt[id];
        if (evt.ec == EC_TERM) continue;
        bool positive = (evtid2base_attempt_rate[id] > 0);
        bool has_legs = !(evtid2outlegs[id].empty() && evtid2inlegs[id].empty());
        for (auto* ls : { &evtid2outlegs[id], &evtid2inlegs[id] }) for (auto& el : *ls) {
            rate a = inflt_attempt_rate[(size_t) id * n_at_slots + el.slot];
            if (!(a < INFINITY)) throw "legs cannot have infinite attempt rates for establishment or act events";
            if (a > 0) positive = true;
        }
        evtid2summary[id] = positive ? 1 : 0;
        if (!(positive && has_legs)) continue;
        evtid_is_leg_summary[id] = 1;
        auto& lsd = evtid2lsd[id];
        vector<int> u0(evtid2outlegs[id].size(), 0), v0(evtid2inlegs[id].size(), 0);
        lsd.H1[u0] = et2n[evt.et1];
        lsd.H3[v0] = et2n[evt.et3];
        rate A; probunits B;
        lsd.T = (double) et2n[evt.et1] * et2n[evt.et3] * class_rate(id, u0, v0, A, B);
        for (auto* ls : { &evtid2outlegs[id], &evtid2inlegs[id] }) for (auto& el : *ls) {
            legkey2summaries[el.legkey].push_back(id);
            if (inflt_attempt_rate[(size_t) id * n_at_slots + el.slot] > 0) lt2sampler[channel_link_type(evt, legkeys[el.legkey])];
        }
    }
}

/** Remember the summary event of a leg summary event type (called by init_events()). */
void init_leg_summary (int evt_id, event& summary_ev)
{
    evtid2summary_ev[evt_id] = summary_ev;
}

/** \returns the total attempt rate of the proposal channels of a summary event with leg influences. */
double leg_summary_channel_total (int evt_id)
{
    auto& evt = evtid2evt[evt_id];
    double N1 = et2n[evt.et1], N3 = et2n[evt.et3];
    double total = evtid2base_attempt_rate[evt_id] * N1 * N3;
    for (auto& el : evtid2outlegs[evt_id]) {
        rate a = inflt_attempt_rate[(size_t) evt_id * n_at_slots + el.slot];
        if (a > 0) total += a * lt2sampler.at(channel_link_type(evt, legkeys[el.legkey])).v.size() * N3;
    }
    for (auto& el : evtid2inlegs[evt_id]) {
        rate a = inflt_attempt_rate[(size_t) evt_id * n_at_slots + el.slot];
        if (a > 0) total += a * lt2sampler.at(channel_link_type(evt, legkeys[el.legkey])).v.size() * N1;
    }
    return total;
}

/** Draw the pair of an attempt of a summary event with leg influences:
 *  choose a channel proportionally to its attempt rate, then a uniformly random pair (basic channel)
 *  or a uniformly random link of the channel's leg type and a uniformly random entity at the other end.
 */
void draw_leg_summary_pair (int evt_id, entity& e1, entity& e3)
{
    auto& evt = evtid2evt[evt_id];
    double N1 = et2n[evt.et1], N3 = et2n[evt.et3];
    double r = uniform(random_variable) * leg_summary_channel_total(evt_id);
    r -= evtid2base_attempt_rate[evt_id] * N1 * N3;
    if (r >= 0) {
        for (int side = 0; side < 2; side++) {
            for (auto& el : (side == 0) ? evtid2outlegs[evt_id] : evtid2inlegs[evt_id]) {
                rate a = inflt_attempt_rate[(size_t) evt_id * n_at_slots + el.slot];
                if (!(a > 0)) continue;
                auto& smp = lt2sampler.at(channel_link_type(evt, legkeys[el.legkey]));
                double w = a * smp.v.size() * ((side == 0) ? N3 : N1);
                if ((r < w) && !smp.v.empty()) {
                    auto& l = smp.v[min((size_t) (uniform(random_variable) * smp.v.size()), smp.v.size() - 1)];
                    if (side == 0) { e1 = l.first; e3 = random_entity(evt.et3); }
                    else { e1 = random_entity(evt.et1); e3 = l.second; }
                    return;
                }
                r -= w;
            }
        }
    }
    // basic channel (also as a fallback against rounding):
    e1 = random_entity(evt.et1);
    e3 = random_entity(evt.et3);
}

/** \returns the class rate sum over all pairs of a leg summary event type, computed from the histograms. */
double leg_summary_T (int evt_id)
{
    return evtid2lsd[evt_id].T;
}

/** \returns the class rate sum over all pairs of a leg summary event type, computed from the histograms. */
double recompute_leg_summary_T (int evt_id)
{
    auto& lsd = evtid2lsd[evt_id];
    double T = 0.0;
    rate A; probunits B;
    for (auto& [u, n1] : lsd.H1) for (auto& [v, n3] : lsd.H3) T += (double) n1 * n3 * class_rate(evt_id, u, v, A, B);
    return T;
}

/** Process the change of a leg degree of an entity. */
static void leg_change (entity e, int k, int delta, const tricllink& l)
{
    const auto& lk = legkeys[k];
    // the degree vectors of e before the change are needed below, so collect the affected summary event types first:
    e2legdeg[(size_t) e * legkeys.size() + k] += delta;

    // individually scheduled events of e (except those about the changing link itself):
    int slot = lk.out ? at_slot(lk.rat, lk.et2, NO_RAT) : at_slot(NO_RAT, lk.et2, lk.rat);
    static vector<event> evs;
    evs.clear();
    for (auto& ev : (lk.out ? e2src_evs[e] : e2tgt_evs[e])) evs.push_back(ev);
    for (auto& ev : evs) {
        if ((ev.e1 == l.e1) && (ev.rat13 == l.rat13) && (ev.e3 == l.e3)) continue;
        int evt_id = evt_id_of(ev.ec, e2et[ev.e1], ev.rat13, e2et[ev.e3]);
        if (evt_id < 0) continue;
        size_t idx = (size_t) evt_id * n_at_slots + slot;
        int j = evtid_at2infl[idx];
        if (j < 0) continue;
        auto evd_ = &ev2data.at(ev);
        unschedule_event(ev, evd_, evt_id);
        if (delta > 0) {
            evd_->n_infl[j] += delta;
            for (int i = 0; i < delta; i++) {
                add_attempt_contribution(evd_, inflt_attempt_rate[idx]);
                add_probunits_contribution(evd_, inflt_delta_probunits[idx]);
            }
        } else {
            assert ((int) evd_->n_infl[j] >= -delta);
            evd_->n_infl[j] += delta;
            for (int i = 0; i < -delta; i++) {
                remove_attempt_contribution(evd_, inflt_attempt_rate[idx]);
                remove_probunits_contribution(evd_, inflt_delta_probunits[idx]);
            }
        }
        _schedule_event(ev, evd_, evt_id);
    }

    // summary event types whose covered pairs depend on this leg degree:
    for (int id : legkey2summaries[k]) {
        auto& evt = evtid2evt[id];
        if (e2et[e] != (lk.out ? evt.et1 : evt.et3)) continue;
        auto& lsd = evtid2lsd[id];
        const auto& ls = lk.out ? evtid2outlegs[id] : evtid2inlegs[id];
        vector<int> w_new = lk.out ? outdegs(id, e) : indegs(id, e), w_old = w_new;
        for (size_t i = 0; i < ls.size(); i++) if (ls[i].legkey == k) w_old[i] -= delta;
        rate A; probunits B;
        // change of the sum over all pairs:
        double dT = 0.0;
        for (auto& [w, n] : (lk.out ? lsd.H3 : lsd.H1)) {
            double lam_new = lk.out ? class_rate(id, w_new, w, A, B) : class_rate(id, w, w_new, A, B);
            double lam_old = lk.out ? class_rate(id, w_old, w, A, B) : class_rate(id, w, w_old, A, B);
            dT += n * (lam_new - lam_old);
            if (compute_gradient) {
                if (lk.out) { register_class_gradient(id, w_new, w, n); register_class_gradient(id, w_old, w, -n); }
                else { register_class_gradient(id, w, w_new, n); register_class_gradient(id, w, w_old, -n); }
            }
        }
        auto& H = lk.out ? lsd.H1 : lsd.H3;
        if (--H[w_old] == 0) H.erase(w_old);
        H[w_new]++;
        lsd.T += dT;
        // pairs of e that are not covered: their shares were subtracted with the old degree vector of e:
        auto correct = [&](entity other) {
            vector<int> w = lk.out ? indegs(id, other) : outdegs(id, other);
            if (other == e) w = lk.out ? indegs(id, e) : outdegs(id, e);
            double lam_new = lk.out ? class_rate(id, w_new, w, A, B) : class_rate(id, w, w_new, A, B);
            double lam_old = lk.out ? class_rate(id, w_old, w, A, B) : class_rate(id, w, w_old, A, B);
            total_finite_effective_rate -= lam_new - lam_old;
            if (compute_gradient) {
                if (lk.out) { register_class_gradient(id, w_new, w, -1); register_class_gradient(id, w_old, w, +1); }
                else { register_class_gradient(id, w, w_new, -1); register_class_gradient(id, w, w_old, +1); }
            }
        };
        if (evt.et1 == evt.et3) correct(e);  // the pair of e with itself
        if (evt.ec == EC_EST) {  // existing links
            if (lk.out) {
                for (auto& ol : e2outs[e]) if ((ol.rat_out == evt.rat13) && (e2et[ol.e_target] == evt.et3)) correct(ol.e_target);
            } else {
                for (auto& il : e2ins[e]) if ((il.rat_in == evt.rat13) && (e2et[il.e_source] == evt.et1)) correct(il.e_source);
            }
        }
        for (auto& ev : (lk.out ? e2src_evs[e] : e2tgt_evs[e])) {  // individually scheduled events
            if ((ev.ec == evt.ec) && (ev.rat13 == evt.rat13) && (e2et[ev.e1] == evt.et1) && (e2et[ev.e3] == evt.et3)) correct(lk.out ? ev.e3 : ev.e1);
        }
        if (!evtid_dirty[id]) { evtid_dirty[id] = 1; dirty_list.push_back(id); }
    }
}

/** Update all leg data after a link was added (delta = +1, after inserting it into the leg sets)
 *  or deleted (delta = -1, after removing it).
 */
void on_link_change (const tricllink& l, int delta)
{
    if (!track_legs || (l.rat13 == RT_ID)) return;
    auto et1 = e2et[l.e1], et3 = e2et[l.e3];
    // link samplers for the leg channels:
    auto sit = lt2sampler.find({ et1, l.rat13, et3 });
    if (sit != lt2sampler.end()) {
        auto& smp = sit->second;
        if (delta > 0) {
            smp.pos[l] = smp.v.size();
            smp.v.push_back({ l.e1, l.e3 });
        } else {
            int i = smp.pos.at(l), last = (int) smp.v.size() - 1;
            if (i != last) { smp.v[i] = smp.v[last]; smp.pos[{ smp.v[i].first, l.rat13, smp.v[i].second }] = i; }
            smp.v.pop_back();
            smp.pos.erase(l);
        }
        for (int id = 0; id < n_evt_ids; id++) {
            if (!evtid_is_leg_summary[id] || evtid_dirty[id]) continue;
            for (int side = 0; side < 2; side++) for (auto& el : (side == 0) ? evtid2outlegs[id] : evtid2inlegs[id]) {
                if (channel_link_type(evtid2evt[id], legkeys[el.legkey]) == link_type{ et1, l.rat13, et3 } && !evtid_dirty[id]) {
                    evtid_dirty[id] = 1; dirty_list.push_back(id);
                }
            }
        }
    }
    // outgoing leg of the source entity and incoming leg of the target entity:
    int k = slot2legkey[at_slot(l.rat13, et3, NO_RAT)];
    if (k >= 0) leg_change(l.e1, k, delta, l);
    k = slot2legkey[at_slot(NO_RAT, et1, l.rat13)];
    if (k >= 0) leg_change(l.e3, k, delta, l);
    // reschedule the summary events whose total rate or channels changed:
    for (int id : dirty_list) {
        evtid_dirty[id] = 0;
        auto& sev = evtid2summary_ev[id];
        auto it = ev2data.find(sev);
        if (it == ev2data.end()) continue;  // (not yet created during initialization)
        unschedule_event(sev, &it->second, id);
        it->second.attempt_rate = leg_summary_channel_total(id);
        _schedule_event(sev, &it->second, id);
    }
    dirty_list.clear();
}
