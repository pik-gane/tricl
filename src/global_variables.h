/*
 * global_variables.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_GLOBAL_VARIABLES_H
#define INC_GLOBAL_VARIABLES_H

// during debugging, you may sometimes want to set the following to true:
#define COUNT_ALL_ANGLES false

#include "data_model.h"


// PARAMETERS

extern string config_yaml_filename;
extern bool verbose, quiet, debug;
extern string diagram_fileprefix;
extern timepoint max_t;
extern long int max_n_events;
extern unsigned seed;
extern unordered_map<relationship_or_action_type, string> gexf_filename;

// structure parameters:

extern unordered_map<entity_type, label> et2label;
extern unordered_map<string, entity_type> label2et;
extern unordered_map<entity_type, entity> et2n; // no. of entities by type
extern int n_rats; // number of distinct relationship_or_action_types
extern unordered_map<relationship_or_action_type, label> rat2label; // relationship_or_action_type labels are verbs or math symbols
extern unordered_map<string, relationship_or_action_type> label2rat;
extern unordered_map<relationship_or_action_type, bool> r_is_action_type; // whether it is an action type
extern unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv; // inverse relations
// any preregistered entities:
extern unordered_map<entity, label> e2label;
extern unordered_map<string, entity> label2e;
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

extern unordered_map<event_type, rate> evt2base_attempt_rate; // basic attempt rate
extern unordered_map<influence_type, rate> inflt2attempt_rate; // note that legs may only attempt removals, not additions!
extern unordered_map<event_type, double> evt2left_tail, evt2right_tail; // tail indices for sigmoid function transforming probunits into probabilities
extern unordered_map<event_type, probunits> evt2base_probunits; // basic success probability unit
extern unordered_map<influence_type, probunits> inflt2delta_probunits; // change in success probunit (legs can also influence additions)

// gexf parameters:
extern unordered_map<entity_type, double> et2gexf_size, et2gexf_a;
extern unordered_map<entity_type, string> et2gexf_shape;
extern unordered_map<entity_type, int> et2gexf_r, et2gexf_g, et2gexf_b;
extern unordered_map<relationship_or_action_type, double> rat2gexf_thickness, rat2gexf_a;
extern unordered_map<relationship_or_action_type, string> rat2gexf_shape;
extern unordered_map<relationship_or_action_type, int> rat2gexf_r, rat2gexf_g, rat2gexf_b;

// DATA:

// arrays of parameters:
extern rate _inflt2attempt_rate[MAX_N_INFLT];
extern probunits _inflt2delta_probunit[MAX_N_INFLT];
extern entity_type e2et[MAX_N_E];

// derived constants:
extern entity max_e; // largest entity id in use
extern set<entity> es;
extern unordered_map<entity_type_pair, set<relationship_or_action_type>> ets2relations;  // possible relations
extern unordered_map<event, rate> ev2max_sp; // max possible success probability of summary events

// variable data:

extern timepoint current_t;
extern long int n_events;
extern event current_ev;
extern event_data* current_evd_;
// network state:
extern unordered_map<entity_type, vector<entity>> et2es;  // kept to equal inverse of v2vt
extern unordered_map<entity, outleg_set> e2outs;
extern unordered_map<entity, inleg_set> e2ins;
extern unordered_map<link_type, long int> lt2n;
extern long int n_links; // total no. of current (non-id.) links incl. inverse relationships
extern long int n_angles; // total no. of angles (excl. RT_ID)

// event data:
extern unordered_map<event, event_data> ev2data;
extern map<timepoint, event> t2be;  // kept to equal inverse of ev2data.t

#endif
