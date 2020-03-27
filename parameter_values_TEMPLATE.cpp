/*
 * parameter_values_TEMPLATE.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 *
 */

/*
 * USAGE: copy to some new file  parameter_values_XXX.cpp,  edit, and insert filename into CMakeList.txt
 *
 */

/* DEFINITIONS:
 * entity type: a type et of concrete or abstract entity (see below), e.g. "an individual", "a news channel", "a social group", "an opinion", or "an infection state"
 * entity: any concrete or abstract object e that can stand in a relationship with other entities, e.g. "John", "the BBC", "catholics", "Elvis lives", or "infected with Dengue"
 * relationship type: any concrete or abstract type rt of directed relationship r two entities can stand in, e.g. "is friends with", "is a subscriber of", "belongs to", "holds", or "is".
 *                    the special relationship type RT_ID encodes the identity of an entity with itself.
 * action type: any type of thing at that can happen at a singular time point between two entities, e.g. "kisses" or "utters that"
 * act: a pair of entities e1,e3 plus an action type at13 plus a time-point t.
 *      we say the act "occurs" at time t.
 * link type: a pair of entity types plus a relationship or action type, e.g. "an individual - is a subscriber of - a news channel" or "an individual - utters - an opinion"
 * link: a pair of entities plus a relationship or action type, e.g. "John - is a subscriber of - the BBC" .
 *       if it has a relationship type, a link "exists" as long as the entities stand in the respective relationship.
 *       if it has an action type, a link "flashes" whenever a corresponding act occurs.
 *       the source entity is called "e1", the destination entity "e3" (because of angles, see below)
 * relationship: a link of relationship type
 * action: a link of action type (not to be confused with an act!)
 * impact of an action: a nonnegative real number that is increased by one whenever the link flashes and decays exponentially at a certain rate.
 *       a link's impact determines the link's influence on adjacent events.
 * basic event: an event class plus a link, e.g.
 *  - "establishment of the link 'John - is a subscriber of - the BBC'" (event class EC_EST)
 *  - "termination of the link 'John - is friends with - Alice'" (event class EC_TERM)
 *  - "occurrence of the act 'Alice - utters that - Elvis lives'" (event class EC_ACT)
 * leg: a relationship or action type rat plus an entity e. a leg can influence an adjacent termination event
 * angle: a pair of relationship or action types rat12,rat23 plus a "middle" entity e2. an angle can influence any adjacent event.
 *        an angle can also encode a leg if either rat12 or rat23 equals NO_RT
 * influence: a basic event plus a leg or angle that influences it (or NO_ANGLE if the basic event happens spontaneously)
 * attempt rate of an influence: the probability rate with which the influence attempts the corresponding event
 * success_probunit of a basic event: probunit of the probability that an attempted event actually happens.
 */

/* ABBREVIATIONS USED IN NAMING:
 * x2y: map mapping xs to ys
 * e[t]: entity (aka node) [type]
 * ev[t]: event [type]
 * infl[t]: influence [type]
 * l[t]: link [type]
 * rat: relationship or action type
 * t: timepoint
 * trailing _: pointer
 * trailing _it: pointer used as iterator
 */

#include "data_model.h"
#include "global_variables.h"


////////////////
// EXECUTION: //
////////////////

const bool verbose = false, quiet = false, debug = false;

// total model time simulated:
timepoint max_t = 1.0;

// where to write network evolution to:
string gexf_filename = "/tmp/tricl.gexf";


////////////////
// STRUCTURE: //
////////////////

// entity type labels:
unordered_map<entity_type, label> et2label = {
//      {et, label}
        {1, "user"},
        {2, "message"},
        {3, "tag"},
        {4, "url"}
};

// no. of entities by type:
unordered_map<entity_type, entity> et2n = {
//      {et, n}
        {1, 1000},
        {2, 100},
        {3, 10},
        {4, 10}
};

// relationship_or_action_type labels (verbs or math symbols):
unordered_map<relationship_or_action_type, label> rat2label = {
//      {rat, label}
        {RT_ID, "="},  // don't change this line (identity relation)
        {2, "–"},
        {3, "→"},
        {4, "posts"},
        {5, "is posted by"}
};

// whether relationship_or_action_type is an action type (true) or a relationship type (false):
unordered_map<relationship_or_action_type, bool> r_is_action_type = {
//      {rat, bool}
        {RT_ID, false},  // don't change this line (identity relation)
        {2, false},
        {3, false},
        {4, true},
        {5, true}
};

// inverse relations:
unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv = {
//      {rat, inv rat}
        {RT_ID, RT_ID},  // don't change this line (identity relation)
        {2, 2}, // 2 is symmetric
        {3, NO_RAT}, // 3 has no inverse
        {4, 5}, // 4's inverse is 5
        {5, 4}, // must be specified though it's redundant
};


////////////////////
// INITIAL STATE: //
////////////////////

// any preregistered entities:
unordered_map<entity, entity_type> e2et = { // their types
//      { e, et }
};
unordered_map<entity, label> e2label = { // their labels
//      { e, label }
};

// any preregistered initial links:
set<link> initial_links = {
//      { e1, rat13, e3 }
};

// random initial links:

// block model:
unordered_map<entity_type, int> et2n_blocks = { // no. of blocks by entity-type (if wanted)
//      { et, n }
        { 1, 30 },
};
unordered_map<link_type, probability> lt2initial_prob_within = { // within-block linking probability by link type
//      { { et1, rat13, et3 }, p }
        { { 1, 2, 1 } , 0.9 },
};
unordered_map<link_type, probability> lt2initial_prob_between = { // between-blocks linking probability by link type
//      { { et1, rat13, et3 }, p }
        { { 1, 2, 1 } , 0.1 },
};

// random geometric model:
unordered_map<entity_type, int> et2dim = { // no. of spatial dimensions by entity-type
//      { et, dim }
        { 2, 2 },
};
unordered_map<link_type, probability> lt2spatial_decay = { // rate of exponential decay of link probability with euclidean distance by link type
//      { { et1, rat13, et3 }, rate }
        { { 2, 2, 2 } , 2.3 },
};


///////////////
// DYNAMICS: //
///////////////

// attempt rates by influence type:
unordered_map<influence_type, rate> inflt2attempt_rate = {
//      { { { ec, et1, rat13, .et3 }, NO_ANGLE }, ar }  // basic attempt rate
//      { { { ec, et1, rat13, .et3 }, { rat12, et2, rat23 } }, ar }  // additional attempt rate due to an angle from e1 to e3
//      { { { EC_TERM, et1, rat13, .et3 }, { rat12, et2, NO_RAT } }, ar }  // termination only: additional attempt rate due to an out-leg at e1
//      { { { EC_TERM, et1, rat13, .et3 }, { NO_RAT, et2, rat23 } }, ar }  // termination only: additional attempt rate due to an in-leg at e3
        { { { EC_EST, 1, 2, 1 }, NO_ANGLE }, 0.001 }, //
        { { { EC_EST, 1, 2, 1 }, { 2, 1, 2 } }, .01 },
        { { { EC_TERM, 1, 2, 1 }, NO_ANGLE }, 1.0 },
};

// basic success probability units:
unordered_map<event_type, probunit> evt2base_probunit = {
//      { { ec, et1, rat13, et3 }, p },
        { { EC_EST, 1, 2, 1 }, INFINITY },
        { { EC_TERM, 1, 2, 1 }, 0 },
};

// change in success probability units:
unordered_map<influence_type, probunit> inflt2delta_probunit = {
//      { { { ec, et1, rat13, .et3 }, { rat12, et2, rat23 } }, pu }  // additional probunits due to an angle from e1 to e3
//      { { { ec, et1, rat13, .et3 }, { rat12, et2, NO_RAT } }, pu }  // additional probunits due to an out-leg at e1
//      { { { ec, et1, rat13, .et3 }, { NO_RAT, et2, rat23 } }, pu }  // additional probunits due to an in-leg at e3
        { { { EC_TERM, 1, 2, 1 }, { 2, 1, 2 } }, -0.3 }, // angles hamper termination
        { { { EC_TERM, 1, 2, 1 }, { 2, 1, NO_RAT } }, +20 }, // legs boost deletion
        { { { EC_TERM, 1, 2, 1 }, { NO_RAT, 1, 2 } }, +20 }, // legs boost deletion
        { { { EC_EST, 1, 2, 1 }, { 2, 1, NO_RAT } }, -10.0 }, // legs hamper establishment
};

// END


/*
 * YAML:
 *  entity-types:
 *   ids:    [1,     2]
 *   labels: [state, individual]
 *   ns:     [3,     10]
 *  relationship types:
 *   ids:      [1,  2]
 *   labels:   [is, meets]
 *   inverses: [-,  2]
 *  action types:
 *   ids: ...
 *   labels: ...
 *  entities:
 *   ids:    [0,           1,        2]
 *   types:  [1,           1,        1]
 *   labels: [susceptible, infected, recovered]
 *  attempts:
 *   "+ 2 -2> 2": 0.1
 *   "- 2 -2> 2": 0.01
 *   "+ 2 -2> 2 -2> 2 -2>": 1.0
 *  base_probs:
 *   "+ 2 -2> 2": 0.5
 *   "- 2 -2> 2": 0.5
 *  influences:
 *   (as in attempts)
 */
