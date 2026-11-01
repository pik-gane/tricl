/** Actions with decaying impacts and events with time-varying rates (see actions.h).
 *
 *  \file
 */

#include <queue>

#include "global_variables.h"
#include "probability.h"
#include "event.h"
#include "actions.h"

bool any_actions = false;
unordered_map<event, tv_data> ev2tv;
vector<vector<double>> evtid2infl_decay;
vector<char> rat_is_action;

/** Impact of an action at a reference time. */
struct impact_data
{
    double I_ref;      ///< Impact at t_ref
    timepoint t_ref;   ///< Reference time (time of the last act or of the creation)
};
static unordered_map<tricllink, impact_data> l2impact;

/** Scheduled expiry of an action (valid only if the action's reference time still equals t_ref). */
struct expiry
{
    timepoint t;
    entity e1;
    relationship_or_action_type rat13;
    entity e3;
    timepoint t_ref;
    bool operator> (const expiry& o) const {
        if (t != o.t) return t > o.t;
        if (e1 != o.e1) return e1 > o.e1;
        if (rat13 != o.rat13) return rat13 > o.rat13;
        return e3 > o.e3;
    }
};
static std::priority_queue<expiry, vector<expiry>, std::greater<expiry>> expiries;

// 8-point Gauss-Legendre nodes and weights on [-1, 1]:
static const double GL_X[8] = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363 };
static const double GL_W[8] = {  0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
                                 0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763 };

double rat_decay (relationship_or_action_type rat)
{
    if ((rat == NO_RAT) || (rat >= rat_is_action.size()) || !rat_is_action[rat]) return 0.0;
    return rat2decay.at(rat);
}

double link_weight (entity e1, relationship_or_action_type rat, entity e3)
{
    if ((rat >= rat_is_action.size()) || !rat_is_action[rat]) return 1.0;
    auto it = l2impact.find({ e1, rat, e3 });
    if (it == l2impact.end()) return 0.0;
    return it->second.I_ref * exp(- rat2decay.at(rat) * (current_t - it->second.t_ref));
}

/** Set up the decay rates of the influences (after init_types()). */
void init_actions ()
{
    rat_is_action.assign(n_rat_slots, 0);
    for (auto& [rat, a] : r_is_action_type) if (a && (rat < (relationship_or_action_type) n_rat_slots)) rat_is_action[rat] = 1;
    any_actions = !initial_impacts.empty();
    evtid2infl_decay.assign(n_evt_ids, {});
    for (int id = 0; id < n_evt_ids; id++) {
        if (rat_is_action[evtid2evt[id].rat13]) any_actions = true;
        for (size_t j = 0; j < evtid2infl_slots[id].size(); j++) {
            int slot = evtid2infl_slots[id][j];
            auto rat12 = (relationship_or_action_type) (slot / (n_et_slots * n_rat_slots));
            auto rat23 = (relationship_or_action_type) (slot % n_rat_slots);
            double d = rat_decay(rat12) + rat_decay(rat23);
            evtid2infl_decay[id].push_back(d);
            if (d > 0) {
                any_actions = true;
                size_t idx = (size_t) id * n_at_slots + slot;
                if (!std::isfinite(inflt_attempt_rate[idx]) || !std::isfinite(inflt_delta_probunits[idx]))
                    throw "angles involving actions must have finite attempt rates and probability units";
            }
        }
    }
}

/** \returns the decaying part of an event's influences and its weights at the current time. */
static inline void decayed_weights (const tv_data& tv, int evt_id, double* W_now)
{
    double dt = current_t - tv.t_ref;
    const auto& dec = evtid2infl_decay[evt_id];
    for (size_t j = 0; j < dec.size(); j++) W_now[j] = (dec[j] > 0) ? tv.W[j] * exp(- dec[j] * dt) : 0.0;
}

rate tv_rate (const event& ev, const event_data* evd_, int evt_id, rate& A, probunits& B, double* W_now)
{
    double W[MAX_INFL_PER_EVT];
    double* w = W_now ? W_now : W;
    decayed_weights(ev2tv.at(ev), evt_id, w);
    A = total_attempt_rate(evd_);
    B = total_success_probunits(evd_);
    const auto& slots = evtid2infl_slots[evt_id];
    const auto& dec = evtid2infl_decay[evt_id];
    for (size_t j = 0; j < dec.size(); j++) if ((dec[j] > 0) && (w[j] > 0)) {
        size_t idx = (size_t) evt_id * n_at_slots + slots[j];
        A += inflt_attempt_rate[idx] * w[j];
        if (std::isfinite(B)) B += inflt_delta_probunits[idx] * w[j];
    }
    if ((A <= 0) || (B == -INFINITY)) return 0.0;
    return A * probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
}

rate tv_bound (const event& ev, const event_data* evd_, int evt_id, rate& A_bound, probunits& B_bound)
{
    double w[MAX_INFL_PER_EVT];
    decayed_weights(ev2tv.at(ev), evt_id, w);
    A_bound = total_attempt_rate(evd_);
    B_bound = total_success_probunits(evd_);
    const auto& slots = evtid2infl_slots[evt_id];
    const auto& dec = evtid2infl_decay[evt_id];
    for (size_t j = 0; j < dec.size(); j++) if ((dec[j] > 0) && (w[j] > 0)) {
        size_t idx = (size_t) evt_id * n_at_slots + slots[j];
        A_bound += inflt_attempt_rate[idx] * w[j];  // (attempt rates are non-negative, so A decreases over time)
        if (std::isfinite(B_bound)) B_bound += max(0.0, inflt_delta_probunits[idx] * w[j]);
    }
    if ((A_bound <= 0) || (B_bound == -INFINITY)) return 0.0;
    return A_bound * probunits2probability(B_bound, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
}

/** Evaluate an event's rate (and the integrand of the gradient) at offset s after the reference time. */
static inline double integrand (const tv_data& tv, const event_data* evd_, int evt_id, double s, double* grad)
{
    const auto& slots = evtid2infl_slots[evt_id];
    const auto& dec = evtid2infl_decay[evt_id];
    double w[MAX_INFL_PER_EVT];
    rate A = total_attempt_rate(evd_);
    probunits B = total_success_probunits(evd_);
    for (size_t j = 0; j < dec.size(); j++) {
        w[j] = (dec[j] > 0) ? tv.W[j] * exp(- dec[j] * s) : 0.0;
        if (w[j] > 0) {
            size_t idx = (size_t) evt_id * n_at_slots + slots[j];
            A += inflt_attempt_rate[idx] * w[j];
            if (std::isfinite(B)) B += inflt_delta_probunits[idx] * w[j];
        }
    }
    if (B == -INFINITY) return 0.0;
    double sigma = probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
    if (grad) {
        double dsigma = probunits2probability_derivative(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
        int k = evtid2base_attempt_param[evt_id];
        if (k >= 0) grad[k] += sigma;
        k = evtid2base_probunits_param[evt_id];
        if (k >= 0) grad[k] += A * dsigma;
        const auto& ia = evtid2infl_attempt_param[evt_id];
        const auto& ib = evtid2infl_probunits_param[evt_id];
        for (size_t j = 0; j < dec.size(); j++) {
            double n = (dec[j] > 0) ? w[j] : (double) evd_->n_infl[j];
            if (n == 0) continue;
            if (ia[j] >= 0) grad[ia[j]] += n * sigma;
            if (ib[j] >= 0) grad[ib[j]] += n * A * dsigma;
        }
    }
    return A * sigma;
}

/** Accumulate the integral of an event's rate (and of its gradient) since its reference time into the log-likelihood,
 *  and move its reference time to the current time.
 */
void settle (const event& ev, event_data* evd_, int evt_id)
{
    auto it = ev2tv.find(ev);
    if (it == ev2tv.end()) return;
    auto& tv = it->second;
    double L = current_t - tv.t_ref;
    if (!(L > 0)) { tv.t_ref = current_t; return; }
    const auto& dec = evtid2infl_decay[evt_id];
    if (evd_->tv) {
        double dmax = 0.0;
        for (size_t j = 0; j < dec.size(); j++) if (tv.W[j] > 0) dmax = max(dmax, dec[j]);
        // pieces of length at most 1/dmax (on which 8-point Gauss-Legendre quadrature is accurate to about 1e-13);
        // after 40 such pieces the decaying parts are below e^-40 of their initial values and the rest is constant:
        double h = (dmax > 0) ? min(L, 1.0 / dmax) : L;
        static vector<double> grad;
        double* g = nullptr;
        if (compute_gradient) { grad.assign(params.size(), 0.0); g = grad.data(); }
        double integral = 0.0, s = 0.0;
        for (int piece = 0; (piece < 40) && (s < L); piece++) {
            double len = min(h, L - s);
            for (int i = 0; i < 8; i++) {
                double x = s + len * (GL_X[i] + 1) / 2, wq = GL_W[i] * len / 2;
                if (g) {
                    static vector<double> gi;
                    gi.assign(params.size(), 0.0);
                    integral += wq * integrand(tv, evd_, evt_id, x, gi.data());
                    for (size_t k = 0; k < gi.size(); k++) g[k] += wq * gi[k];
                } else {
                    integral += wq * integrand(tv, evd_, evt_id, x, nullptr);
                }
            }
            s += len;
        }
        if (s < L) {  // constant remainder
            double len = L - s;
            if (g) {
                static vector<double> gi;
                gi.assign(params.size(), 0.0);
                integral += len * integrand(tv, evd_, evt_id, L, gi.data());
                for (size_t k = 0; k < gi.size(); k++) g[k] += len * gi[k];
            } else {
                integral += len * integrand(tv, evd_, evt_id, L, nullptr);
            }
        }
        cumulative_logl -= integral;
        if (g) for (size_t k = 0; k < grad_exposure.size(); k++) grad_exposure[k] += g[k];
    }
    for (size_t j = 0; j < dec.size(); j++) if (dec[j] > 0) tv.W[j] *= exp(- dec[j] * L);
    tv.t_ref = current_t;
}

/** Settle all events with time-varying rates (at the end of a simulation or replay). */
void settle_all ()
{
    for (auto& [ev, tv] : ev2tv) {
        auto it = ev2data.find(ev);
        if (it != ev2data.end()) settle(ev, &it->second, evt_id_of(ev.ec, e2et[ev.e1], ev.rat13, e2et[ev.e3]));
    }
}

/** Add the gradient of the log rate of a time-varying event that happens now to \ref grad_event_terms. */
void add_tv_event_log_gradient (const event& ev, const event_data* evd_, int evt_id)
{
    if (!compute_gradient) return;
    double w[MAX_INFL_PER_EVT];
    rate A; probunits B;
    tv_rate(ev, evd_, evt_id, A, B, w);
    if (!(A > 0) || (B == -INFINITY)) return;
    double sigma = probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]),
           dsigma = probunits2probability_derivative(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
    double dls = (sigma > 0) ? dsigma / sigma : 0.0;
    int k = evtid2base_attempt_param[evt_id];
    if (k >= 0) grad_event_terms[k] += 1 / A;
    k = evtid2base_probunits_param[evt_id];
    if (k >= 0) grad_event_terms[k] += dls;
    const auto& dec = evtid2infl_decay[evt_id];
    const auto& ia = evtid2infl_attempt_param[evt_id];
    const auto& ib = evtid2infl_probunits_param[evt_id];
    for (size_t j = 0; j < dec.size(); j++) {
        double n = (dec[j] > 0) ? w[j] : (double) evd_->n_infl[j];
        if (n == 0) continue;
        if (ia[j] >= 0) grad_event_terms[ia[j]] += n / A;
        if (ib[j] >= 0) grad_event_terms[ib[j]] += n * dls;
    }
}

timepoint next_expiry ()
{
    while (!expiries.empty()) {
        auto& x = expiries.top();
        auto it = l2impact.find({ x.e1, x.rat13, x.e3 });
        if ((it != l2impact.end()) && (it->second.t_ref == x.t_ref)) return x.t;
        expiries.pop();  // outdated entry
    }
    return INFINITY;
}

/** Remove an action whose impact has fallen to its threshold, and update the adjacent events. */
static void remove_action (const tricllink& l)
{
    auto it = l2impact.find(l);
    if (it == l2impact.end()) return;
    double I = it->second.I_ref * exp(- rat2decay.at(l.rat13) * (current_t - it->second.t_ref));
    l2impact.erase(it);
    e2outs[l.e1].erase({ .rat_out = l.rat13, .e_target = l.e3 });
    e2ins[l.e3].erase({ .e_source = l.e1, .rat_in = l.rat13 });
    lt2n[{ e2et[l.e1], l.rat13, e2et[l.e3] }]--;
    n_links--;
    update_adjacent_angles(l.e1, l.rat13, l.e3, -1, -I);
}

void process_expiries_until (timepoint t)
{
    while (true) {
        timepoint tx = next_expiry();
        if (!(tx <= t)) break;
        auto x = expiries.top();
        expiries.pop();
        advance_time(max(current_t, tx));
        if (debug) cout << "at t=" << current_t << " action " << e2label[x.e1] << " " << rat2label[x.rat13] << " " << e2label[x.e3] << " expires" << endl;
        remove_action({ x.e1, x.rat13, x.e3 });
        auto inv = rat2inv.at(x.rat13);
        if (inv != NO_RAT) remove_action({ x.e3, inv, x.e1 });  // (the mirrored action has the same impact)
    }
}

/** Increase the impact of one action by an amount (creating the action if necessary). */
static void jump (const tricllink& l, double amount)
{
    double delta = rat2decay.at(l.rat13), theta = rat2threshold.at(l.rat13);
    auto it = l2impact.find(l);
    double I;
    int dcount;
    if (it != l2impact.end()) {
        I = it->second.I_ref * exp(- delta * (current_t - it->second.t_ref)) + amount;
        it->second = { I, current_t };
        dcount = 0;
    } else {
        I = amount;
        if (!(I >= theta)) return;
        l2impact[l] = { I, current_t };
        e2outs[l.e1].insert({ .rat_out = l.rat13, .e_target = l.e3 });
        e2ins[l.e3].insert({ .e_source = l.e1, .rat_in = l.rat13 });
        lt2n[{ e2et[l.e1], l.rat13, e2et[l.e3] }]++;
        n_links++;
        dcount = 1;
    }
    expiries.push({ current_t + log(I / theta) / delta, l.e1, l.rat13, l.e3, current_t });
    update_adjacent_angles(l.e1, l.rat13, l.e3, dcount, amount);
}

void act_on_link (const tricllink& l, double amount)
{
    jump(l, amount);
    auto inv = rat2inv.at(l.rat13);
    if (inv != NO_RAT) jump({ l.e3, inv, l.e1 }, amount);
}

/** Create the named initial actions with their impacts (during initialization). */
void add_initial_actions ()
{
    for (auto& [l, I] : initial_impacts) {
        if (l2impact.count(l) > 0) continue;  // (already created as the mirror of another initial action)
        act_on_link(l, I);
    }
}
