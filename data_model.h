/*
 * data_model.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_DATA_MODEL_H
#define INC_DATA_MODEL_H

#include <math.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>

using namespace std;

// TYPES AND STRUCTS:

// basic types:

typedef size_t relationship_or_action_type;
#define NO_RAT 0 // missing value
#define RT_ID 1 // identity relationship

typedef int entity; // actual entities >= 0. can be used to store entity types as negative numbers.
#define E(e) e
#define MAX_N_E 1000

#define EPS 1e-20;

typedef double timepoint;
typedef double rate;
typedef double probunit;
typedef double probability;

enum event_class { EC_EST, EC_TERM, EC_ACT };
extern unordered_map<event_class, string> ec2label;

typedef string label;

// complex types:

typedef entity entity_type; // >=1 so that -entity_type can be stored in entity fields

// bit sizes for hashs:
#define E_BITS 12
#define RAT_BITS 2
#define ET_BITS 2
// must have 2*E_BITS+RAT_BITS <= 30 and RAT_BITS+ET_BITS <= 10

struct entity_type_pair
{
    const entity_type et1;  // type of source entity e1
    const entity_type et3;  // type of destination entity e3

    friend bool operator==(const entity_type_pair& left, const entity_type_pair& right) {
        return (left.et1 == right.et1
                && left.et3 == right.et3);
    }
};
template <> struct std::hash<entity_type_pair> {
    size_t operator()(const entity_type_pair& etp) const {
        return (etp.et1 ^ (etp.et3 << ET_BITS));
    }
};

struct link_type
{
    const entity_type et1;  // type of source entity e1
    const relationship_or_action_type rat13;  // relation of link e1 -> e3
    const entity_type et3;  // type of destination entity e3

    friend bool operator==(const link_type& left, const link_type& right) {
        return (left.et1 == right.et1
                && left.rat13 == right.rat13
                && left.et3 == right.et3);
    }
    friend ostream& operator<<(ostream& os, const link_type& lt);
};
template <> struct std::hash<link_type> {
    size_t operator()(const link_type& lt) const {
        return (lt.et1 ^ (lt.rat13 << ET_BITS) ^ (lt.et3 << (ET_BITS+RAT_BITS)));
    }
};

struct event_type // TODO: inherit from link_type?
{
    const event_class ec;
    const entity_type et1;  // type of source entity e1
    const relationship_or_action_type rat13;  // relation of link e1 -> e3
    const entity_type et3;  // type of destination entity e3

    friend bool operator==(const event_type& left, const event_type& right) {
        return (left.ec == right.ec
                && left.et1 == right.et1
                && left.rat13 == right.rat13
                && left.et3 == right.et3);
    }
    friend ostream& operator<<(ostream& os, const event_type& e);
};
template <> struct std::hash<event_type> {
    size_t operator()(const event_type& evt) const {
        return (evt.ec ^ (evt.et1 << 2) ^ (evt.rat13 << (2+ET_BITS)) ^ (evt.et3 << (2+ET_BITS+RAT_BITS)));
    }
};

struct angle_type
{
    const relationship_or_action_type rat12;  // relationship_or_action_type of link e1 -> e2, or NO_RT if in-leg
    const entity_type et2;  // type of middle entity e3
    const relationship_or_action_type rat23;  // relationship_or_action_type of link e2 -> e3, or NO_RT if out-leg

    friend bool operator==(const angle_type& left, const angle_type& right) {
        return (left.rat12 == right.rat12
                && left.et2 == right.et2
                && left.rat23 == right.rat23);
    }
};
template <> struct std::hash<angle_type> {
    size_t operator()(const angle_type& at) const {
        return (at.rat12 ^ (at.et2 << RAT_BITS) ^ (at.rat23 << (RAT_BITS+ET_BITS)));
    }
};

const angle_type NO_ANGLE = { .rat12 = NO_RAT, .et2 = 0, .rat23 = NO_RAT }; // signifying spontaneous events

struct influence_type
{
    const event_type evt;
    const angle_type at; // or NO_ANGLE if event is spontaneous

    friend bool operator==(const influence_type& left, const influence_type& right) {
        return (left.evt == right.evt
                && left.at == right.at);
    }
};
#define INFLT(inflt) ((size_t)inflt.evt.ec ^ ((size_t)inflt.evt.et1 << 2) ^ ((size_t)inflt.evt.rat13 << (2+ET_BITS)) ^ ((size_t)inflt.evt.et3 << (2+ET_BITS+RAT_BITS)) ^ ((size_t)inflt.at.rat12 << (2+2*ET_BITS+RAT_BITS)) ^ ((size_t)inflt.at.et2 << (2+2*ET_BITS+2*RAT_BITS)) ^ ((size_t)inflt.at.rat23 << (2+3*ET_BITS+2*RAT_BITS)))
#define MAX_N_INFLT (1 << (2+3*ET_BITS+3*RAT_BITS))
template <> struct std::hash<influence_type> {
    size_t operator()(const influence_type& inflt) const {
        return (INFLT(inflt));
    }
};

struct leg // occurs in out-link and in-link sets
{
    const entity e;
    const relationship_or_action_type r;

    friend bool operator==(const leg& left, const leg& right) {
        return (left.e == right.e
                && left.r == right.r);
    }
    friend bool operator<(const leg& left, const leg& right) {
        return ((left.e < right.e) || ((left.e == right.e) && (left.r < right.r)));
    }
};
template <> struct std::hash<leg> {
    size_t operator()(const leg& l) const {
        return (l.e ^ (l.r << E_BITS));
    }
};

typedef set<leg> leg_set;

struct link
{
    const entity e1;  // source entity (or, if <0, source entity type of spontaneous addition)
    const relationship_or_action_type rat13;  // relationship (!) type of link e1 -> e3
    const entity e3;  // destination entity (or, if <0, dest. entity type of spontaneous addition)

    friend bool operator==(const link& left, const link& right) {
        return (left.e1 == right.e1
                && left.rat13 == right.rat13
                && left.e3 == right.e3);
    }
};
template <> struct std::hash<link> {
    size_t operator()(const link& l) const {
        return (l.e1 ^ (l.rat13 << E_BITS) ^ (l.e3 << (E_BITS+RAT_BITS)));
    }
};

struct event // TODO: inherit from link?
{
    event_class ec;
    entity e1;  // source entity (or, if <0, source entity type of spontaneous addition)
    relationship_or_action_type rat13;  // relationship_or_action_type of link e1 -> e3
    entity e3;  // destination entity (or, if <0, dest. entity type of spontaneous addition)

    friend bool operator==(const event& left, const event& right) {
        return (left.ec == right.ec
                && left.e1 == right.e1
                && left.rat13 == right.rat13
                && left.e3 == right.e3);
    }
    friend ostream& operator<<(ostream& os, const event& ev);
};
template <> struct std::hash<event> {
    size_t operator()(const event& ev) const {
        return (ev.ec ^ (ev.e1 << 2) ^ (ev.e3 << (2+E_BITS)) ^ (ev.rat13 << (2+E_BITS+E_BITS)));
    }
};

struct event_data
{
    int n_angles;
    rate attempt_rate;
    probunit success_probunit;
    timepoint t = -INFINITY;
};

struct angle
{
    relationship_or_action_type rat12;  // relationship_or_action_type of link e1 -> e2
    entity e2;  // middle entity
    relationship_or_action_type rat23;  // relationship_or_action_type of link e2 -> e3

    friend bool operator==(const angle& left, const angle& right) {
        return (left.rat12 == right.rat12
                && left.e2 == right.e2
                && left.rat23 == right.rat23);
    }
};
template <> struct std::hash<angle> {
    size_t operator()(const angle& a) const {
        return (a.rat12 ^ (a.e2 << RAT_BITS) ^ (a.rat23 << (RAT_BITS+E_BITS)));
    }
};

typedef vector<angle> angles;

#endif
