// make sure this file is only included once:
#ifndef INC_LEGS_H
#define INC_LEGS_H

/** Influences of legs (links of one endpoint of an event to entities of a given type) on events.
 *
 *  \file
 *
 *  A leg influence makes the attempt rate or the success probability units of an event depend on the number of
 *  links of a given relationship type from its source entity to entities of a given type ("outgoing leg",
 *  [~, rel.type, ent.type] in the config file) or from entities of a given type to its target entity
 *  ("incoming leg", [ent.type, rel.type, ~]), not counting the link the event is about.
 *
 *  Individually scheduled events count their legs in event_data::n_infl like angles and are updated whenever a leg
 *  changes (via the per-entity indices of scheduled events). Pairs covered by a summary event (establishment or act
 *  events without influencing angles) have rates that depend on the leg degrees of their two entities; their total
 *  rate is maintained exactly as
 *
 *      sum over all pairs of the class rate lambda(u(e1), v(e3))  -  sum over the pairs not covered of the same,
 *
 *  where u and v are the vectors of outgoing and incoming leg degrees relevant for the event type,
 *  the first sum is kept via histograms of these degree vectors, and the pairs not covered (equal entities,
 *  existing links, individually scheduled events) are enumerated per entity when a degree changes.
 *  Attempts of summary events are proposed through channels (the basic attempt rate for uniformly drawn pairs,
 *  and per leg influence with positive attempt rate a uniformly drawn link of the leg's type, i.e. a source or
 *  target entity drawn proportionally to its degree) and accepted with probability f(B) / p_max,
 *  which yields exactly the rate A * f(B) of every covered pair (thinning).
 */

#include "global_variables.h"
#include "probability.h"

/** Description of a leg type that influences at least one event type. */
struct legkey_info
{
    bool out;                          ///< Whether this is an outgoing leg of the source entity (else an incoming leg of the target entity)
    relationship_or_action_type rat;   ///< Relationship type of the leg
    entity_type et2;                   ///< Type of the other entity of the leg
};

/** A leg influence of an event type. */
struct evt_leg
{
    int j;       ///< Influence index within the event type
    int legkey;  ///< Index in \ref legkeys
    int slot;    ///< Angle type slot (with NO_RAT for the missing relationship type)
};

/** Data of a summary event type whose covered pairs are influenced by legs. */
struct leg_summary_data
{
    map<vector<int>, long> H1;  ///< Number of source entities by vector of outgoing leg degrees
    map<vector<int>, long> H3;  ///< Number of target entities by vector of incoming leg degrees
    double T = 0.0;             ///< Sum of the class rates over all ordered pairs of entities of the two types
};

/** Uniform sampling of the links of a link type (used for the leg channels of summary events). */
struct link_sampler
{
    vector<pair<entity, entity>> v;      ///< Source and target entities of the links
    unordered_map<tricllink, int> pos;   ///< Position of each link in v
};

extern bool track_legs;                                   ///< Whether any leg influences exist (then the leg data below is maintained)
extern vector<legkey_info> legkeys;                       ///< All leg types that influence some event type
extern vector<int> slot2legkey;                           ///< Index in legkeys by angle type slot, or -1
extern vector<int> e2legdeg;                              ///< Leg degree by entity * legkeys.size() + leg type index
extern vector<vector<evt_leg>> evtid2outlegs, evtid2inlegs;  ///< Leg influences by event type id
extern vector<char> evtid_is_leg_summary;                 ///< Whether the summary event of an event type has leg influences
extern vector<leg_summary_data> evtid2lsd;                ///< Leg summary data by event type id
extern unordered_map<link_type, link_sampler> lt2sampler; ///< Link samplers by link type (only for link types used by leg channels)
extern vector<unordered_set<event>> e2src_evs, e2tgt_evs; ///< Individually scheduled events by source and by target entity (only if track_legs)
extern vector<int> evtid2summary;                         ///< Whether an event type has a summary event (1) or not (0), by event type id

/** \returns the leg degree of an entity for a leg type. */
inline int legdeg (entity e, int k)
{
    return e2legdeg[(size_t) e * legkeys.size() + k];
}

/** Compute the attempt rate and success probunits of a pair covered by a summary event from the leg degree vectors
 *  u (outgoing legs of the source) and v (incoming legs of the target), and return its rate.
 */
inline rate class_rate (int evt_id, const vector<int>& u, const vector<int>& v, rate& A, probunits& B)
{
    A = evtid2base_attempt_rate[evt_id];
    B = evtid2base_probunits[evt_id];
    const auto& ol = evtid2outlegs[evt_id];
    for (size_t i = 0; i < ol.size(); i++) if (u[i] > 0) {
        size_t idx = (size_t) evt_id * n_at_slots + ol[i].slot;
        A += u[i] * inflt_attempt_rate[idx];
        B += u[i] * inflt_delta_probunits[idx];
    }
    const auto& il = evtid2inlegs[evt_id];
    for (size_t i = 0; i < il.size(); i++) if (v[i] > 0) {
        size_t idx = (size_t) evt_id * n_at_slots + il[i].slot;
        A += v[i] * inflt_attempt_rate[idx];
        B += v[i] * inflt_delta_probunits[idx];
    }
    if ((A == 0) || (B == -INFINITY)) return 0.0;
    return A * probunits2probability(B, evtid2left_tail[evt_id], evtid2right_tail[evt_id], evtid2scale[evt_id]);
}

/** \returns the vector of outgoing leg degrees of a source entity relevant for an event type. */
inline vector<int> outdegs (int evt_id, entity e)
{
    const auto& ol = evtid2outlegs[evt_id];
    vector<int> u(ol.size());
    for (size_t i = 0; i < ol.size(); i++) u[i] = legdeg(e, ol[i].legkey);
    return u;
}

/** \returns the vector of incoming leg degrees of a target entity relevant for an event type. */
inline vector<int> indegs (int evt_id, entity e)
{
    const auto& il = evtid2inlegs[evt_id];
    vector<int> v(il.size());
    for (size_t i = 0; i < il.size(); i++) v[i] = legdeg(e, il[i].legkey);
    return v;
}

/** \returns the rate of a pair if it is covered by the summary event of its event type (including leg influences). */
inline rate covered_pair_rate (int evt_id, entity e1, entity e3, rate& A, probunits& B)
{
    if (!evtid_is_leg_summary[evt_id]) {
        A = evtid2base_attempt_rate[evt_id];
        B = evtid2base_probunits[evt_id];
        return evtid2summary_single_er[evt_id];
    }
    return class_rate(evt_id, outdegs(evt_id, e1), indegs(evt_id, e3), A, B);
}

void register_class_gradient (int evt_id, const vector<int>& u, const vector<int>& v, double weight);
void add_class_log_gradient (int evt_id, const vector<int>& u, const vector<int>& v);
void add_pair_share (int evt_id, entity e1, entity e3);
void subtract_pair_share (int evt_id, entity e1, entity e3);
void add_covered_event_log_gradient (int evt_id, entity e1, entity e3);
void index_event (const event& ev);
void unindex_event (const event& ev);
void init_legs ();
void init_leg_summary (int evt_id, event& summary_ev);
double leg_summary_channel_total (int evt_id);
void draw_leg_summary_pair (int evt_id, entity& e1, entity& e3);
void on_link_change (const tricllink& l, int delta);
double recompute_leg_summary_T (int evt_id);

#endif
