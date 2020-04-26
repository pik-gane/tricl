// make sure this file is only included once:
#ifndef INC_GLOBAL_VARIABLES_H
#define INC_GLOBAL_VARIABLES_H

/** Definition of global data storage
 *
 * \file
 *
 * See \ref data_model.h for an explanation of the overall data architecture.
 */

#include "data_model.h"

// during debugging, you may sometimes want to set the following to true:
#define COUNT_ALL_ANGLES false

// CONSTANT DATA:

extern unordered_set<entity> es;               ///< Set of all entities
extern entity max_e;                           ///< Largest entity id in use
extern entity_type e2et[MAX_N_E];              ///< Entity type by entity (stored in an array for performance)
extern unordered_map<entity_type, vector<entity>> et2es;  ///< Inverse of e2et
extern unordered_map<entity, label> e2label;   ///< Labels of entities
extern unordered_map<string, entity> label2e;  ///< Inverse map of \ref e2label

// Options and parameters:

extern string config_yaml_filename; ///< Name of (or path to) main config file
extern bool only_output_logl;       ///< Whether to restrict output to log-likelihood information
extern bool debug;                  ///< Whether to output debug messages
extern bool silent;                 ///< Whether to suppress all output except what was requested explicitly
extern bool quiet;                  ///< Whether to suppress most output
extern bool verbose;                ///< Whether to output more detailed information
extern string diagram_fileprefix;   ///< Prefix of name of (or path to) generated diagram files
extern timepoint max_t;             ///< Maximal model time to simulate until
extern long int max_n_events;       ///< Max. no. events to simulate before stopping
extern unsigned seed;               ///< Random seed (if 0, generate a random seed)
extern unordered_map<relationship_or_action_type, string> gexf_filename;  ///< Names of (or paths to) generated gexf (or gexf.gz) files by relationship or action type

// structure parameters:

extern unordered_map<entity_type, label> et2label;                    ///< Entity type labels (typically nouns)
extern unordered_map<string, entity_type> label2et;                   ///< Inverse map of \ref et2label
extern unordered_map<entity_type, entity> et2n;                       ///< No. of entities by type
extern int n_rats;                                                    ///< No. of distinct relationship or action types
extern unordered_map<relationship_or_action_type, label> rat2label;   ///< Relationship or action type labels (typically verbs in the 3rd person singular, or math symbols)
extern unordered_map<string, relationship_or_action_type> label2rat;  ///< Inverse map of \ref label2rat
extern unordered_map<relationship_or_action_type, bool> r_is_action_type; ///< Whether relationship or action type is an action type (not implemented yet)
extern unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv; ///< Inverse type of a relationship or action type (e.g. the inverse of "follows" would be "is followed by", the inverse of "meets" would be "meets"). If NO_RAT, inverse has no individual label
// any preregistered initial links:
extern set<tricllink> initial_links;                                       ///< Set of named initial links
// random initial links:
// block model:
extern unordered_map<entity_type, int> et2n_blocks;                   ///< No. of blocks for random block model by entity type
extern unordered_map<link_type, probability> lt2initial_prob_within;  ///< Within-block link probability for random block model by link type
extern unordered_map<link_type, probability> lt2initial_prob_between; ///< Between-block link probability for random block model by link type
// random geometric model:
extern unordered_map<entity_type, int> et2dim;                        ///< No. of spatial dimensions for random geometric model by entity type
extern unordered_map<link_type, probability> lt2spatial_decay;        ///< Rate of exponential decay of link probability for random geometric model by link type

// dynamic parameters:

extern unordered_set<event_type> possible_evts;                                    ///< Types of events that may occur at all
extern unordered_map<event_type, rate> evt2base_attempt_rate;            ///< Basic attempt rate by event type
extern unordered_map<influence_type, rate> inflt2attempt_rate;           ///< Additional attempt rate by influence type
extern rate _inflt2attempt_rate[MAX_N_INFLT];                            ///< Redundant copy of \ref inflt2attempt_rate as array
extern unordered_map<event_type, double> evt2left_tail,                  ///< Left tail index for sigmoid function probunits2probability(), >=0
                                         evt2right_tail;                 ///< Right tail index for sigmoid function probunits2probability(), >= 0
extern unordered_map<event_type, probunits> evt2base_probunits;          ///< Basic success probability units by event type
extern unordered_map<influence_type, probunits> inflt2delta_probunits;   ///< Change in success probunits by influence type
extern probunits _inflt2delta_probunits[MAX_N_INFLT];                    ///< Redundant copy of \ref inflt2delta_probunits as array
extern unordered_map<entity_type_pair, unordered_set<relationship_or_action_type>> ets2relations;  ///< Possible relationship or action types by entity type pair
extern unordered_map<event, rate> summary_ev2max_success_probability;    ///< Maximal possible success probability of summary events
extern unordered_map<event_type, rate> summary_evt2single_effective_rate;  ///< Effective rate of a single instance of a summary event

// gexf parameters:
extern unordered_map<entity_type, double> et2gexf_size,                         ///< Node size for gexf file by entity type
                                          et2gexf_a;                            ///< Node color alpha value for gexf file by entity type, 0...1
extern unordered_map<entity_type, string> et2gexf_shape;                        ///< Node shape for gexf file by entity type
extern unordered_map<entity_type, int> et2gexf_r,                               ///< Node color red value for gexf file by entity type, 0...255
                                       et2gexf_g,                               ///< Node color green value for gexf file by entity type, 0...255
                                       et2gexf_b;                               ///< Node color blue value for gexf file by entity type, 0...255
extern unordered_map<relationship_or_action_type, double> rat2gexf_thickness,   ///< Edge line thickness for gexf file by relationship or action type
                                                          rat2gexf_a;           ///< Edge color alpha value for gexf file by relationship or action type, 0...1
extern unordered_map<relationship_or_action_type, string> rat2gexf_shape;       ///< Edge shape for gexf file by relationship or action type
extern unordered_map<relationship_or_action_type, int> rat2gexf_r,              ///< Edge color red value for gexf file by relationship or action type, 0...255
                                                       rat2gexf_g,              ///< Edge color green value for gexf file by relationship or action type, 0...255
                                                       rat2gexf_b;              ///< Edge color blue value for gexf file by relationship or action type, 0...255

// VARIABLE DATA:

extern timepoint current_t;       ///< Current model time point
extern timepoint last_dt;         ///< Time between last and current event
extern long int n_events;         ///< No. of events that occurred so far
extern event current_ev;          ///< Current event
extern event_data* current_evd_;  ///< Pointer to event data of current event
extern map<timepoint, event> t2ev;                ///< Current schedule of events, inverse of ev2data[ev].t. CAUTION: this needs to be an ordered container type
extern unordered_map<event, event_data> ev2data;  ///< Data of all currently scheduled events

// log-likelihood computation:
extern double cumulative_logl;    ///< Cumulative log-likelihood of evolution from initial state to current_t
extern rate summary_ev2effective_rate;  ///< Current effective rates of summary events
extern rate total_finite_effective_rate; ///< Current total effective rate of all events (without infinite rates)
extern int n_infinite_effective_rates; ///< No. of currently scheduled events with infinite effective rate
extern event_data sure_evd;

// network state:
extern unordered_map<entity, outleg_set> e2outs;  ///< Set of current outlegs by source entity
extern unordered_map<entity, inleg_set> e2ins;    ///< Set of current inlegs by source entity (redundant, but essential for performance)
extern unordered_map<link_type, long int> lt2n;   ///< No. of current (non-id.) links by type incl. inverse relationships
extern long int n_links;                          ///< Total no. of current (non-id.) links incl. inverse relationships
extern long int n_angles;                         ///< Total no. of current (non-id.) angles that may influence at least one event

#endif
