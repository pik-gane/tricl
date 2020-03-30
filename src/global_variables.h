/*
 * global_variables.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_GLOBAL_VARIABLES_H
#define INC_GLOBAL_VARIABLES_H

#include <string>
#include <map>
#include <unordered_map>

#include "data_model.h"

using namespace std;

// PARAMETERS

extern string config_yaml_filename;

extern bool verbose, quiet, debug;

extern string gexf_filename;

extern timepoint max_t;

// structure parameters:

extern unordered_map<entity_type, label> et2label;
extern unordered_map<entity_type, entity> et2n; // no. of entities by type
extern unordered_map<relationship_or_action_type, label> rat2label; // relationship_or_action_type labels are verbs or math symbols
extern unordered_map<relationship_or_action_type, bool> r_is_action_type; // whether it is an action type
extern unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv; // inverse relations
// any preregistered entities:
extern unordered_map<entity, entity_type> e2et;
extern unordered_map<entity, label> e2label;
// any preregistered initial links:
extern set<link> initial_links;
// random initial links:
// block model:
extern unordered_map<entity_type, int> et2n_blocks;
extern unordered_map<link_type, probability> lt2initial_prob_within;
extern unordered_map<link_type, probability> lt2initial_prob_between;
// random geometric model:
extern unordered_map<entity_type, int> et2dim;
extern unordered_map<link_type, probability> lt2spatial_decay;

// dynamic parameters:

extern unordered_map<influence_type, rate> inflt2attempt_rate; // note that legs may only attempt removals, not additions!
extern unordered_map<event_type, probunit> evt2base_probunit; // basic success probability unit
extern unordered_map<influence_type, probunit> inflt2delta_probunit; // change in success probunit (legs can also influence additions)


// DATA:

// arrays of parameters:
extern rate _inflt2attempt_rate[MAX_N_INFLT];
extern probunit _inflt2delta_probunit[MAX_N_INFLT];
extern entity_type _e2et[MAX_N_E];

// derived constants:
extern entity max_e; // largest entity id in use
extern set<entity> es;
extern unordered_map<entity_type_pair, set<relationship_or_action_type>> ets2relations;  // possible relations

// variable data:

extern timepoint current_t;
extern event current_ev;
extern event_data* current_evd_;
// network state:
extern unordered_map<entity_type, vector<entity>> et2es;  // kept to equal inverse of v2vt
extern unordered_map<entity, leg_set> e2outs;
extern unordered_map<entity, leg_set> e2ins;
extern int n_links; // total no. of current (non-id.) links incl. inverse relationships
// event data:
extern unordered_map<event, event_data> ev2data;
extern map<timepoint, event> t2be;  // kept to equal inverse of ev2data.t

#endif
