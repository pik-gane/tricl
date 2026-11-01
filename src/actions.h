// make sure this file is only included once:
#ifndef INC_ACTIONS_H
#define INC_ACTIONS_H

/** Actions: links of an action type, whose impact jumps by 1 at every act and decays exponentially.
 *
 *  \file
 *
 *  An act is an event of class EC_ACT on a pair of entities and an action type. It does not change whether a link
 *  exists but increases the impact I of the action (the link of that action type between the two entities, and of
 *  its mirrored link for symmetric or inverse action types) by 1. Between acts, I decays as I(t) = I(t0) e^(-delta (t - t0))
 *  with the decay rate delta of the action type, and the action ceases to exist when I falls below the threshold
 *  theta of the action type (at a deterministic time). This threshold is part of the model, so that the set of
 *  existing links stays finite and the simulation and likelihood remain exact.
 *
 *  An angle (or the identity angle of the action itself, e.g. [~, posts, message, =, ~]) whose legs include actions
 *  contributes to the influence sums of adjacent events with the product of the impacts of its legs (relationship
 *  links count 1), so that its weight decays at the sum of the decay rates of its action legs. Events with such
 *  "decaying" influences have time-varying rates
 *
 *      A(t) = A_static + sum_j a_j W_j(t),   B(t) = B_static + sum_j b_j W_j(t),   W_j(t) = W_j(t0) e^(-delta_j (t - t0)),
 *
 *  where the static parts contain the basic values, legs and angles of relationship links only.
 *  Since attempt rates are non-negative, A(t) <= A(t0) and B(t) <= B_static + sum_j max(0, b_j W_j(t0)) between
 *  updates, so the rate tree holds the bound A(t0) f(B_bound), and an event drawn from the tree happens with
 *  probability (rate at the current time) / bound (thinning; otherwise its bound is refreshed).
 *  The integral of the rate (and of its gradient) over time, needed for the log-likelihood, is accumulated
 *  lazily per event by Gauss-Legendre quadrature whenever the event is updated ("settled").
 */

#include "global_variables.h"

extern bool any_actions;                               ///< Whether the model has action types with dynamics or initial actions
extern unordered_map<event, tv_data> ev2tv;            ///< Time-varying data of the events with decaying influences
extern vector<vector<double>> evtid2infl_decay;        ///< Decay rate of each influence of each event type (0 for static influences)
extern vector<char> rat_is_action;                     ///< Whether a relationship or action type is an action type, by id

/** \returns the decay rate of a relationship or action type (0 for relationship types and NO_RAT). */
double rat_decay (relationship_or_action_type rat);

/** \returns the current weight of a link in angles: its impact for actions, 1 otherwise. */
double link_weight (entity e1, relationship_or_action_type rat, entity e3);

void init_actions ();
void settle (const event& ev, event_data* evd_, int evt_id);
void settle_all ();
rate tv_rate (const event& ev, const event_data* evd_, int evt_id, rate& A, probunits& B, double* W_now);
rate tv_bound (const event& ev, const event_data* evd_, int evt_id, rate& A_bound, probunits& B_bound);
void add_tv_event_log_gradient (const event& ev, const event_data* evd_, int evt_id);
timepoint next_expiry ();
void process_expiries_until (timepoint t);
void act_on_link (const tricllink& l, double amount);
void process_one_expiry ();
extern unordered_map<relationship_or_action_type, vector<pair<relationship_or_action_type, event_class>>> rat2triggered;  ///< Action types whose acts are triggered by establishment or termination events of a relationship type
void perform_triggered_acts (event_class ec, entity e1, relationship_or_action_type rat13, entity e3);
void add_initial_actions ();

#endif
