// make sure this file is only included once:
#ifndef INC_DATA_MODEL_H
#define INC_DATA_MODEL_H

/** The main data model
 *
 * \file
 *
 * Data architecture
 * -------------------
 *
 * Most data is kept in unordered maps (what would be dictionaries in python)
 * whose keys are entity types, entities, link types, links, legs, angles,
 * event types, events, influence types and influences.
 * the latter types of things are encoded as structs,
 * most of which are basically tuples of integer ids.
 *
 * For fast access to map entries, these are automatically turned into integer
 * hashs via the hash structs at the end of this file.
 *
 * Some data (which is accessed most often) is instead kept in vectors
 * whose indices are either entities (which are ints)
 * or integer hash values of influence types (constructed via the macro INFLT).
 *
 * All hashs are constructed as logical ORs of properly shifted ids,
 * hence valid ids are restricted by the respective numbers of bits reserved
 * for this id in the hash.
 *
 * This hash construction is governed by the bit size macros
 * #E_BITS, #ET_BITS, and #RAT_BITS,
 * which could be adapted in dependence on system architecture,
 * but must fulfil the following constraints:
 *   2 + 2 * #E_BITS + #RAT_BITS <= no. of bits in size_t (32 or 64)
 *   2^#E_BITS + 2^(6 + 3 * #RAT_BITS + 3 * #ET_BITS) <= available memory bytes
 */

// the following choices seem adequate for 64 bit system and >= 1 GB available memory:
#define E_BITS 20   ///< No. of bits used for entities --> max. 1 mio. entities
#define ET_BITS 4   ///< No. of bits used for entity types --> max. 16 entity types
#define RAT_BITS 4  ///< No. of bits used for relationship or action types --> max. 16 relationship or action types

#define MAX_N_E ((1<<E_BITS)-1)  ///< resulting max. no. of entities


// commonly used includes are all included here:

#include <assert.h>  // used in DEBUG mode only

// io:
#include <fstream>
#include <iostream>

// datatypes:
#include <string>
#include <vector>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>

// other:
#include <math.h>


/// We use our own namespace to avoid clashes with 3rdparty names, e.g. "link"
namespace tricl {

// make common members of std available without needing to specify "std::" each time:
using std::ofstream;
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;
using std::vector;
using std::unordered_set;
using std::set;
using std::unordered_map;
using std::map;
using std::max;
using std::min;
using std::round;



// TYPES AND STRUCTS AND CORRESPONDING CONSTANTS:


// basic types:

typedef double timepoint;    ///< A point in continuous model time, 0...inf
typedef double probunits;    ///< -inf...inf, will be mapped to probabilities by means of the function \ref probunits2probability()
typedef double probability;  ///< 0...1
typedef double rate;         ///< Probability per time, 0...inf. (inf is used for things happening "immediately")
typedef string label;        ///< A label of something

/** Entities are used to encode all physical and abstract objects of the modeled social dynamics
 *  that can stand in some form of relation to another or can perform some kinds of actions on or with another.
 *
 *  Actual entities have ids >= 1.
 *  In summary events, the datatype "entity" is also used to store entity types as negative numbers.
 */
typedef int entity;

/** Entity types are used to encode the different kinds of entities in a model.
 *
 *  Examples of entity type labels: "user", "message", "opinion", "epidemic status", "religious group", "household", "news channel"
 */
typedef short unsigned int entity_type; ///< >=1 so that -entity_type can be stored in entity fields


/** Relationship types are used to encode the kinds of relationships entities may stand in;
 *  action types are used to encode the kinds of actions entities may perform on or with another (not implemented yet);
 *  they are encoded in a common datatype and distinguished by means of the map rat_is_action_type (not implemented yet).
 *
 *  Both can either be symmetric (undirected) or nonsymmetric (directed).
 *
 *  Examples of labels for symmetric relationship types: "is friends with", "is inconsistent with", "is a neighbour of".
 *
 *  Examples of labels for nonsymmetric relationship types: "follows", "has read", "is in status", "belongs to".
 *
 *  Examples of labels for symmetric action types: "meets", "talks to".
 *
 *  Examples of labels for nonsymmetric action types: "broadcasts", "reads".
 */
typedef size_t relationship_or_action_type;  ///< >= 1

#define NO_RAT 0  ///< Missing value for relationship or action type, used in angles to encode legs
#define RT_ID 1   ///< Relationship type for identity relationship "=", always present


// complex types:

/** Pair of entity types (rarely used)
 */
struct entity_type_pair
{
    const entity_type et1;  ///< Type of source entity e1
    const entity_type et3;  ///< Type of target entity e3

    friend bool operator== (const entity_type_pair& left, const entity_type_pair& right) {
        return (left.et1 == right.et1
                && left.et3 == right.et3);
    }
};

/** A link encodes either the existence of a certain relationship between two entities
 *  (or the cumulative impact of all past actions of a certain type between two entities -- not yet implemented).
 */
struct link
{
    const entity e1;                          ///< Source entity (or, if <0, source entity type of a summary event)
    const relationship_or_action_type rat13;  ///< Type of relationship or action represented by this link
    const entity e3;                          ///< Target entity (or, if <0, target entity type of a summary event)

    friend bool operator== (const link& left, const link& right) {
        return ((left.e1 == right.e1)
                && (left.rat13 == right.rat13)
                && (left.e3 == right.e3));
    }
    // since links are used in sets, they must provide a comparison operator:
    friend bool operator< (const link& left, const link& right) {
        return ((left.e1 < right.e1)
                || ((left.e1 == right.e1) && (left.rat13 < right.rat13))
                || ((left.e1 == right.e1) && (left.rat13 == right.rat13) && (left.e3 < right.e3)));
    }
};

/** The type of a \ref link is given by two entity types and a relationship or action type.
 */
struct link_type
{
    const entity_type et1;                    ///< Type of source entity e1 occurring in this type of link
    const relationship_or_action_type rat13;  ///< Type of relationship or action represented by this type of link
    const entity_type et3;                    ///< Type of target entity e3 occurring in this type of link

    friend bool operator== (const link_type& left, const link_type& right) {
        return (left.et1 == right.et1
                && left.rat13 == right.rat13
                && left.et3 == right.et3);
    }
};

/** (don't confuse with event_type!) The event class encodes what kind of change to the system state an \ref event represents.
 */
enum event_class {
    EC_EST,   ///< Establishment of a relationship
    EC_TERM,  ///< Termination of a relationship
    EC_ACT    ///< Occurrence of an action (not yet implemented)
};

extern unordered_map<event_class, label> ec2label;

/** An event encodes a possible atomic change in the system state involving exactly one link.
 *
 *  For simplicity, the link is not encoded by an attribute of type \ref link
 *  but rather by a copy of its three attributes e1, rat13, e3.
 */
struct event
{
    // in contrast to the other structs, here the members cannot be const since events are used as values of maps:
    event_class ec;                     ///< What kind of change to system state this event represents
    entity e1;                          ///< Source entity (or, if <0, source entity type of a summary event)
    relationship_or_action_type rat13;  ///< Type of relationship or action affected by the event
    entity e3;                          ///< Target entity (or, if <0, target entity type of a summary event)

    friend bool operator== (const event& left, const event& right) {
        return (left.ec == right.ec
                && left.e1 == right.e1
                && left.rat13 == right.rat13
                && left.e3 == right.e3);
    }
};

/** The type of an \ref event is given by an event class, two entity types and a relationship or action type.
 */
struct event_type
{
    const event_class ec;                     ///< What kind of change to system state this type of event represents
    const entity_type et1;                    ///< Type of source entity e1 affected by this type of event
    const relationship_or_action_type rat13;  ///< Type of relationship or action affected by this type of event
    const entity_type et3;                    ///< Type of target entity e3 affected by this type of event

    friend bool operator== (const event_type& left, const event_type& right) {
        return (left.ec == right.ec
                && left.et1 == right.et1
                && left.rat13 == right.rat13
                && left.et3 == right.et3);
    }
};

/** For performance reasons, the mutable data of an \ref event is stored in a separate struct.
 *
 *  These structs appear as values in a map whose key is the corresponding event.
 */
struct event_data
{
    int n_angles = 0;             ///< Current no. of angles influencing this event
    rate attempt_rate;            ///< Current attempt rate of this event
    probunits success_probunits;  ///< Current success probunits of this event
    timepoint t = -INFINITY;      ///< When this event would next happen if the system state does not chance in between
};

/** An inleg represents a leg "incoming" to a target entity.
 *
 *  When an inleg is used, the target entity is clear from the context,
 *  hence it is not named itself in the inleg but only the source entity
 *  and relationship or action type are stored.
 *
 *  Inlegs may influence the attempt rate or success probability of "adjacent" events
 *  (events with the same target entity).
 */
struct inleg
{
    const entity e_source;                     ///< Source entity of the leg. Plays a similar role as \c e2 in an angle.
    const relationship_or_action_type rat_in;  ///< Plays a similar role as \c rat23 in an angle.

    friend bool operator== (const inleg& left, const inleg& right) {
        return (left.e_source == right.e_source
                && left.rat_in == right.rat_in);
    }
    // since inlegs are used in sets, they must provide a comparison operator:
    friend bool operator< (const inleg& left, const inleg& right) {
        // CAUTION: order must be lexicographic with e_other leading (otherwise function leg_intersection() won't work!):
        return ((left.e_source < right.e_source)
                || ((left.e_source == right.e_source) && (left.rat_in < right.rat_in)));
    }
};

/** An outleg represents a leg "outgoing" from a source entity.
 *
 *  When an outleg is used, the source entity is clear from the context,
 *  hence it is not named itself in the outleg but only the target entity
 *  and relationship or action type are stored.
 *
 *  Outlegs may influence the attempt rate or success probability of "adjacent" events
 *  (events with the same source entity).
 */
struct outleg
{
    const relationship_or_action_type rat_out;  ///< Plays a similar role as \c rat12 in an angle.
    const entity e_target;                      ///< Target entity of the leg. Plays a similar role as \c e2 in an angle.

    friend bool operator== (const outleg& left, const outleg& right) {
        return (left.e_target == right.e_target
                && left.rat_out == right.rat_out);
    }
    // since outlegs are used in sets, they must provide a comparison operator:
    friend bool operator< (const outleg& left, const outleg& right) {
        // order must be lexicographic with e_other leading (otherwise leg_intersection won't work!):
        return ((left.e_target < right.e_target) || ((left.e_target == right.e_target) && (left.rat_out < right.rat_out)));
    }
};

// CAUTION: the following must be (ordered) set (NOT unordered_set) for leg_intersection to work!!
typedef set<inleg> inleg_set;    ///< used to store all incoming legs of an entity
typedef set<outleg> outleg_set;  ///< used to store all outgoing legs of an entity

/** An angle represents an indirect connection between a source entity
 *  and a target entity via two links through some "middle" entity.
 *
 *  When an angle is used, the source and target entities are clear from the context,
 *  hence they are not named themselves in the angle but only the middle entity
 *  and relationship or action type are stored.
 *
 *  Angles may influence the attempt rate or success probability of "adjacent" events
 *  (events with the same source and target entities).
 *
 *  The \c angle datatype may also be used to store an \ref inleg by setting \c rat12 to \ref NO_RAT
 *  or an \ref outleg by setting \c rat23 to \ref NO_RAT.
 */
struct angle
{
    // as in events, here the members cannot be const since angles are used as values of maps:
    relationship_or_action_type rat12;  ///< Relationship or action type of link from source (usually \c e1) to middle entity \c e2, or \ref NO_RAT if the angle encodes an \ref inleg
    entity e2;                          ///< Middle entity
    relationship_or_action_type rat23;  ///< Relationship or action type of link from middle entity \c e2 to target (usually \c e3), or \ref NO_RAT if the angle encodes an \ref outleg

    friend bool operator== (const angle& left, const angle& right) {
        return (left.rat12 == right.rat12
                && left.e2 == right.e2
                && left.rat23 == right.rat23);
    }
};

/** The type of an angle is given by the involved entity and relationship or action types
 *  (or \ref NO_RAT if the angle type is actually a leg type).
 */
struct angle_type
{
    const relationship_or_action_type rat12;  ///< Relationship or action type of link from source (usually \c e1) to middle entity \c e2, or \ref NO_RAT if the angle encodes an \ref inleg
    const entity_type et2;                    ///< Type of middle entity
    const relationship_or_action_type rat23;  ///< Relationship or action type of link from middle entity \c e2 to target (usually \c e3), or \ref NO_RAT if the angle encodes an \ref outleg

    friend bool operator== (const angle_type& left, const angle_type& right) {
        return (left.rat12 == right.rat12
                && left.et2 == right.et2
                && left.rat23 == right.rat23);
    }
};

typedef vector<angle> angles;

/** Influence types are used to encode the main dynamical rules of a model.
 *
 *  The attempt rates and success probability units of an event can be influenced by
 *  adjacent angles and legs (not yet implemented), and the amount of influence
 *  is a parameter that depends on the respective entity and relationship or action types.
 *
 *  These parameters are stored in maps that use the influence_type as key.
 */
struct influence_type
{
    const event_type evt;  ///< The type of event influenced by this type of influence
    const angle_type at;   ///< The type of angle (or leg via NO_RAT) influencing that type of event

    friend bool operator== (const influence_type& left, const influence_type& right) {
        return (left.evt == right.evt
                && left.at == right.at);
    }
};

} // end of namespace tricl


using namespace tricl;

// define hashing methods (must be outside namespace tricl):

/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<entity_type_pair> {
    inline size_t operator()(const entity_type_pair& etp) const {
        return (etp.et1 ^ (etp.et3 << ET_BITS));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<link_type> {
    inline size_t operator()(const link_type& lt) const {
        return (lt.et1 ^ (lt.rat13 << ET_BITS) ^ (lt.et3 << (ET_BITS+RAT_BITS)));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<event_type> {
    inline size_t operator()(const event_type& evt) const {
        return (evt.ec ^ (evt.et1 << 2) ^ (evt.rat13 << (2+ET_BITS)) ^ (evt.et3 << (2+ET_BITS+RAT_BITS)));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<inleg> {
    inline size_t operator()(const inleg& l) const {
        return (l.e_source ^ (l.rat_in << E_BITS));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<outleg> {
    inline size_t operator()(const outleg& l) const {
        return (l.e_target ^ (l.rat_out << E_BITS));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<angle> {
    inline size_t operator()(const angle& a) const {
        return (a.rat12 ^ (a.e2 << RAT_BITS) ^ (a.rat23 << (RAT_BITS+E_BITS)));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<angle_type> {
    inline size_t operator()(const angle_type& at) const {
        return (at.rat12 ^ (at.et2 << RAT_BITS) ^ (at.rat23 << (RAT_BITS+ET_BITS)));
    }
};

#define INFLT(inflt) ((size_t)inflt.evt.ec ^ ((size_t)inflt.evt.et1 << 2) ^ ((size_t)inflt.evt.rat13 << (2+ET_BITS)) ^ ((size_t)inflt.evt.et3 << (2+ET_BITS+RAT_BITS)) ^ ((size_t)inflt.at.rat12 << (2+2*ET_BITS+RAT_BITS)) ^ ((size_t)inflt.at.et2 << (2+2*ET_BITS+2*RAT_BITS)) ^ ((size_t)inflt.at.rat23 << (2+3*ET_BITS+2*RAT_BITS)))
#define MAX_N_INFLT (1 << (2+3*ET_BITS+3*RAT_BITS))
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<influence_type> {
    inline size_t operator()(const influence_type& inflt) const {
        return (INFLT(inflt));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<link> {
    inline size_t operator()(const link& l) const {
        return (l.e1 ^ (l.rat13 << E_BITS) ^ (l.e3 << (E_BITS+RAT_BITS)));
    }
};
/** Construct an integer hash for use in maps and sets by adding bit-shifted atteributes:
 */
template <> struct std::hash<event> {
    inline size_t operator()(const event& ev) const {
        return (ev.ec ^ (ev.e1 << 2) ^ (ev.e3 << (2+E_BITS)) ^ (ev.rat13 << (2+E_BITS+E_BITS)));
    }
};

#endif
