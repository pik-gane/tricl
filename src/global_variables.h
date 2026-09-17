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
extern string events_out_filename;  ///< Name of (or path to) csv file to write all performed events to (if "", none is written)
extern string events_in_filename;   ///< Name of (or path to) csv file with events to replay instead of simulating (if "", simulate)
extern bool output_summary;         ///< Whether to output a one-line JSON summary at the end
extern bool compute_gradient;       ///< Whether to compute the gradient of the log-likelihood w.r.t. the model parameters
extern bool dump_parameters;        ///< Whether to only output the model parameters as JSON and exit
extern bool dump_model;             ///< Whether to only output the model structure (types, event types, initial link counts) as JSON after initialization and exit
extern bool scheduling_enabled;     ///< Whether tentative event times are drawn and kept in the schedule (true when simulating, false when replaying)
extern string stats_out_filename;   ///< Name of (or path to) csv file to write link counts by link type to at regular model time intervals (if "", none is written)
extern double stats_every;          ///< Model time interval between rows of the stats file
extern string links_out_filename;   ///< Name of (or path to) csv file to write the time intervals of all links to (if "", none is written)
extern string entities_out_filename; ///< Name of (or path to) csv file to write all entities with their types to (if "", none is written)
extern timepoint max_t;             ///< Maximal model time to simulate until (may be infinite if max_n_events is finite)
extern timepoint never_t;           ///< Finite time point at or after which "never" happening events are formally scheduled (= max_t if finite)
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
extern unordered_set<relationship_or_action_type> inverse_only_rats; ///< Relationship types that were only declared as the inverse of another type (their links are implied by those of the other type and hence not output)
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
extern unordered_map<event_type, rate> evt2base_attempt_rate;            ///< Basic attempt rate by event type (as read from the config file)
extern unordered_map<influence_type, rate> inflt2attempt_rate;           ///< Additional attempt rate by influence type (as read from the config file)
extern unordered_map<event_type, double> evt2left_tail,                  ///< Left tail index for sigmoid function probunits2probability(), >=0
                                         evt2right_tail;                 ///< Right tail index for sigmoid function probunits2probability(), >= 0
extern unordered_map<event_type, probunits> evt2base_probunits;          ///< Basic success probability units by event type (as read from the config file)
extern unordered_map<influence_type, probunits> inflt2delta_probunits;   ///< Change in success probunits by influence type (as read from the config file)
extern bool any_leg_influences;                                          ///< Whether any influence type is a leg (rather than angle) influence with nonzero effect (if not, loops over legs can be skipped)
extern unordered_map<entity_type_pair, unordered_set<relationship_or_action_type>> ets2relations;  ///< Possible relationship or action types by entity type pair

// dense indexing of types for fast parameter lookup in the simulation's hot loops (set up in init_types()):

extern int n_et_slots;   ///< No. of entity type slots (= largest entity type id + 1)
extern int n_rat_slots;  ///< No. of relationship or action type slots (= largest relationship or action type id + 1)
extern int n_at_slots;   ///< No. of angle type slots (= n_rat_slots * n_et_slots * n_rat_slots)
extern int n_evt_ids;    ///< No. of possible event types, which get dense ids 0 ... n_evt_ids-1
extern vector<int> evt_slot2id;                      ///< Dense id of a possible event type by its slot (see evt_slot()), or -1 if the event type cannot occur
extern vector<event_type> evtid2evt;                 ///< Event type by dense id
extern vector<rate> evtid2base_attempt_rate;         ///< Base attempt rate by event type id
extern vector<probunits> evtid2base_probunits;       ///< Base success probunits by event type id
extern vector<double> evtid2left_tail,               ///< Left tail index of the sigmoid by event type id
                      evtid2right_tail,              ///< Right tail index of the sigmoid by event type id
                      evtid2scale;                   ///< Precomputed scale parameter of the sigmoid (see tail2scale()) by event type id
extern vector<rate> evtid2summary_single_er;         ///< Effective rate of a single pair covered by the summary event of this event type (0 if there is no summary event)
extern vector<probability> evtid2summary_max_success_probability;  ///< Upper bound of the success probability used for scheduling the summary event of this event type
extern vector<rate> inflt_attempt_rate;              ///< Additional attempt rate by influence, indexed via inflt_index()
extern vector<probunits> inflt_delta_probunits;      ///< Change in success probunits by influence, indexed via inflt_index()
extern vector<vector<relationship_or_action_type>> ets2rats;  ///< Possible relationship or action types by et1 * n_et_slots + et3

/** \returns the slot of an angle type in the influence tables. */
inline int at_slot (relationship_or_action_type rat12, entity_type et2, relationship_or_action_type rat23)
{
    return ((int) rat12 * n_et_slots + (int) et2) * n_rat_slots + (int) rat23;
}
/** \returns the slot of an event type in \ref evt_slot2id. */
inline int evt_slot (event_class ec, entity_type et1, relationship_or_action_type rat13, entity_type et3)
{
    return (((int) ec * n_et_slots + (int) et1) * n_rat_slots + (int) rat13) * n_et_slots + (int) et3;
}
/** \returns the dense id of an event type, or -1 if the event type cannot occur. */
inline int evt_id_of (event_class ec, entity_type et1, relationship_or_action_type rat13, entity_type et3)
{
    return evt_slot2id[evt_slot(ec, et1, rat13, et3)];
}
/** \returns the index of an influence (event type id plus angle type) in the influence tables. */
inline int inflt_index (int evt_id, relationship_or_action_type rat12, entity_type et2, relationship_or_action_type rat23)
{
    return evt_id * n_at_slots + at_slot(rat12, et2, rat23);
}
/** \returns the effective rate of a single pair covered by the summary event of an event type (0 if there is none). */
inline rate summary_single_er_of (event_class ec, entity_type et1, relationship_or_action_type rat13, entity_type et3)
{
    int evt_id = evt_id_of(ec, et1, rat13, et3);
    return (evt_id >= 0) ? evtid2summary_single_er[evt_id] : 0.0;
}

// influences by event type (an "influence" is an angle type that has a nonzero effect on an event type):

extern vector<vector<int>> evtid2infl_slots;  ///< Angle type slots of the influences of each event type, by event type id (at most MAX_INFL_PER_EVT each)
extern vector<signed char> evtid_at2infl;     ///< Influence index (0 ... MAX_INFL_PER_EVT-1) by event type id * n_at_slots + angle type slot, or -1 if the angle type has no effect on the event type

/** \returns the influence index of an angle type for an event type, or -1 if the angle type has no effect on the event type. */
inline int infl_index_of (int evt_id, relationship_or_action_type rat12, entity_type et2, relationship_or_action_type rat23)
{
    return evtid_at2infl[inflt_index(evt_id, rat12, et2, rat23)];
}

// model parameters, for gradients of the log-likelihood:

/** What kind of model parameter a \ref model_param is. */
enum param_kind {
    PK_BASE_ATTEMPT,     ///< base attempt rate of an event type
    PK_INFL_ATTEMPT,     ///< additional attempt rate of an event type due to an influence
    PK_BASE_PROBUNITS,   ///< base success probunits of an event type
    PK_INFL_PROBUNITS    ///< change of the success probunits of an event type due to an influence
};

/** A model parameter as specified in the dynamics section of the config file. */
struct model_param
{
    param_kind kind;   ///< What kind of parameter this is
    int evt_id;        ///< Event type id the parameter belongs to
    int infl_index;    ///< Influence index within the event type (for PK_INFL_*), or -1
    double value;      ///< Current value
    string label;      ///< Human-readable label, used in JSON output
};

extern vector<model_param> params;                        ///< All model parameters (in the order of their indices)
extern vector<int> evtid2base_attempt_param;              ///< Index in params of the base attempt rate of an event type, or -1 if not specified in the config
extern vector<int> evtid2base_probunits_param;            ///< Index in params of the base probunits of an event type, or -1 if not specified in the config
extern vector<vector<int>> evtid2infl_attempt_param;      ///< Index in params of the attempt rate influence, by event type id and influence index, or -1
extern vector<vector<int>> evtid2infl_probunits_param;    ///< Index in params of the probunits influence, by event type id and influence index, or -1
extern vector<double> grad_event_terms;  ///< Sum over performed events of the gradient of log(event rate), by parameter index
extern vector<double> grad_rate;         ///< Gradient of the current total effective rate, by parameter index
extern vector<double> grad_exposure;     ///< Integral of grad_rate over model time so far, by parameter index (the gradient of the log-likelihood is grad_event_terms - grad_exposure)

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
extern vector<outleg_set> e2outs;  ///< Set of current outlegs by source entity (entity ids are dense, 1 ... max_e)
extern vector<inleg_set> e2ins;    ///< Set of current inlegs by target entity (redundant, but essential for performance)
extern unordered_map<link_type, long int> lt2n;   ///< No. of current (non-id.) links by type incl. inverse relationships
extern long int n_links;                          ///< Total no. of current (non-id.) links incl. inverse relationships
extern long int n_angles;                         ///< Total no. of current (non-id.) angles that may influence at least one event

#endif
