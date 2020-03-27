/*
 * parameter_values_30blocks.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 *
 */

/*
 * simple model with
 * 1000 entities of one type,
 * in 30 disjoint blocks fully connected by a unique symmetric relationship type,
 * slow spontaneous linking,
 * medium-speed linking due to angles (triangle closure),
 * and fast spontaneous unlinking hampered by angles (triangle staibilization),
 * leading to metastable clusters which eventually merge.
 */

#include "data_model.h"
#include "global_variables.h"


////////////////
// EXECUTION: //
////////////////

const bool verbose = false, quiet = false, debug = false;

// total model time simulated:
timepoint max_t = 6.5;

// where to write network evolution to:
string gexf_filename = "/tmp/30blocks.gexf";


////////////////
// STRUCTURE: //
////////////////

// entity type labels:
unordered_map<entity_type, label> et2label = {
//      {et, label}
        {1, "node"}
};

// no. of entities by type:
unordered_map<entity_type, entity> et2n = {
//      {et, n}
        {1, 1000},
};

// relationship_or_action_type labels (verbs or math symbols):
unordered_map<relationship_or_action_type, label> rat2label = {
//      {rat, label}
        {RT_ID, "="},
        {2, "–"},
};

// whether relationship_or_action_type is an action type (true) or a relationship type (false):
unordered_map<relationship_or_action_type, bool> r_is_action_type = {
//      {rat, bool}
        {RT_ID, false},
        {2, false},
};

// inverse relations:
unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv = {
//      {rat, inv rat}
        {RT_ID, RT_ID},  // don't change this line (identity relation)
        {2, 2}, // 2 is symmetric
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
        { { 1, 2, 1 } , 1.0 },
};
unordered_map<link_type, probability> lt2initial_prob_between = { // between-blocks linking probability by link type
//      { { et1, rat13, et3 }, p }
        { { 1, 2, 1 } , 0.0 },
};

// random geometric model:
unordered_map<entity_type, int> et2dim = { // no. of spatial dimensions by entity-type
//      { et, dim }
};
unordered_map<link_type, probability> lt2spatial_decay = { // rate of exponential decay of link probability with euclidean distance by link type
//      { { et1, rat13, et3 }, rate }
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
        { { { EC_EST, 1, 2, 1 }, NO_ANGLE }, 0.001 },
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
};

// END



// other, similar examples would be:
// one symmetric type, N 100, T 200, three blocks, .001 spontaneous, 0.01 est, 3.0 -0.3 term -> metastable
// one symmetric type, N 300, T 100, four blocks, .001 spontaneous, 0.01 est, 10.0 -0.3 term -> metastable
// one symmetric type, N 300, T 15, three blocks, .001 spontaneous, 0.01 est, 10.0 -0.3 term -> metastable, but eventually merging

