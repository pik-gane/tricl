/*
 * tricl.cpp
 *
 *  Created on: Oct 18, 2019
 *      Author: Jobst Heitzig, PIK
 */

// compile with g++ -Wall -std=gnu++17 -O3 -o tricl tricl.cpp -lyaml-cpp
// profile with valgrind --tool=callgrind ./tricl; kcachegrind callgrind.out

/* TODO:
 * add event preclasses EC_A_EST, EC_A_ACT for attempts of generic events, then allow legs to influence their success prob.
 * use const args for optimization
 * take into account influence of legs on addition!
 * add acts
 * read parameters from YAML file
 * output clustering coefficient
 * add particle filtering mode
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

/* OTHER NAMING CONVENTIONS:
 * void functions are named by imperatives (e.g., "do_this")
 * non-void functions are named by nouns (e.g., "most_important_thing")
 */

// INCLUDES:

// for release mode, uncomment this:
#define NDEBUG
#include <assert.h>
//

#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
//#include "yaml-cpp/yaml.h"
#include "yaml.h"

// data type includes:
#include <cfloat>
#include <string>
#include <tuple>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <iterator>

using namespace std;


// TYPES AND STRUCTS:

// basic types:

typedef size_t relationship_or_action_type;
#define NO_RAT 0 // missing value
#define RT_ID 1 // identity relationship

typedef int entity; // actual entities >= 0. can be used to store entity types as negative numbers.
#define E(e) e
#define MAX_N_E 1000

typedef double timepoint;
typedef double rate;
typedef double probunit;
typedef double probability;

enum event_class { EC_EST, EC_TERM, EC_ACT };
unordered_map<event_class, string> ec2label = {
    { EC_EST, "establish" },
    { EC_TERM, "terminate" },
    { EC_ACT, "act" }
};

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

double tail2scale(double tail)
{
    return 1 / (1 + 1 / tail) / pow(1 + log(1 + 1 / tail), tail + 1) / 2;
}

probability probunit2probability (probunit pu)
{
//    return 1.0 / (1.0 + exp(-l));
    double  left_tail = 1, // 1/sigma0 from paper. larger means heavier tail
            right_tail = 1,  // 1/sigma1
            scale = tail2scale(left_tail) + tail2scale(right_tail);
    return (1 / pow(1 + log(1 + exp(- pu / scale) / left_tail), left_tail)
            + 1
            - 1 / pow(1 + log(1 + exp(pu / scale) / right_tail), right_tail) ) / 2;
}


/////////////////
// PARAMETERS: //
/////////////////

// TODO: later read them from file!

const bool verbose = false; bool quiet = false, debug = false;
//const bool verbose = false; bool quiet = true, debug = false;
//const bool verbose = true; bool quiet = false, debug = true;
string gexf_filename = "/tmp/tricl.gexf";

timepoint max_t = 7.0; //100.0 ; //1.0;

// structure parameters:

// examples:
// one symmetric type, N 100, T 200, three blocks, .001 spontaneous, 0.01 est, 3.0 -0.3 term -> metastable
// one symmetric type, N 300, T 100, four blocks, .001 spontaneous, 0.01 est, 10.0 -0.3 term -> metastable
// one symmetric type, N 300, T 15, three blocks, .001 spontaneous, 0.01 est, 10.0 -0.3 term -> metastable, but eventuall merging
// one symmetric type, N 1000, T 7, 30 blocks, .001 spontaneous, 0.01 est, 1.0 -0.3 term -> metastable, slowly merging

unordered_map<entity_type, label> et2label = {
        {1, "node"}
};
unordered_map<entity_type, entity> et2n = { // no. of entities by type
        {1, 1000},
};
unordered_map<relationship_or_action_type, label> rat2label = { // relationship_or_action_type labels are verbs or math symbols
        {RT_ID, "="},
//        {2, "→"},
        {2, "–"},
};
unordered_map<relationship_or_action_type, bool> r_is_action_type = { // whether it is an action type
        {RT_ID, false},
        {2, false},
};
unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv = { // inverse relations
        {RT_ID, RT_ID},
//        {2, NO_RAT},
        {2, 2}, // 2 is symmetric
};
// any preregistered entities:
unordered_map<entity, entity_type> e2et = {
};
unordered_map<entity, label> e2label = {
};
// any preregistered initial links:
set<link> initial_links = {
};
// random initial links:
// block model:
unordered_map<entity_type, int> et2n_blocks = {
        { 1, 30 },
};
unordered_map<link_type, probability> lt2initial_prob_within = {
        { { 1, 2, 1 } , 1.0 },
};
unordered_map<link_type, probability> lt2initial_prob_between = {
        { { 1, 2, 1 } , 0.0 },
};
// random geometric model:
unordered_map<entity_type, int> et2dim = {
//        { 1, 2 },
};
unordered_map<link_type, probability> lt2spatial_decay = {
//        { { 1, 2, 1 } , 4.5 },
};


// dynamic parameters:

unordered_map<influence_type, rate> inflt2attempt_rate = { // note that legs may only attempt removals, not additions!
//        { { { ec, et1, rat13, .et3 }, { rat12, et2, rat23 } }, ar },
        { { { EC_EST, 1, 2, 1 }, NO_ANGLE }, 0.001 },
        { { { EC_EST, 1, 2, 1 }, { 2, 1, 2 } }, .01 },
        { { { EC_TERM, 1, 2, 1 }, NO_ANGLE }, 1.0 },
};

unordered_map<event_type, probunit> evt2base_probunit = {  // basic success probability unit
//        { { ec, et1, rat13, et3 }, p },
        { { EC_EST, 1, 2, 1 }, INFINITY },
        { { EC_TERM, 1, 2, 1 }, 0 },
//        { { EC_TERM, 1, 2, 1 }, -10-20*10 }, // so that with 5 nbs. its negligible but with 6 considerable
};

unordered_map<influence_type, probunit> inflt2delta_probunit = { // change in success probunit (legs can also influence additions)
       { { { EC_TERM, 1, 2, 1 }, { 2, 1, 2 } }, -0.3 }, // angles hamper termination
//        { { { EC_TERM, 1, 2, 1 }, { 2, 1, NO_RAT } }, +20 }, // legs boost deletion
//        { { { EC_TERM, 1, 2, 1 }, { NO_RAT, 1, 2 } }, +20 }, // legs boost deletion
//        { { { EC_EST, 1, 2, 1 }, { 2, 1, NO_RAT } }, -10.0 }, // TODO: legs hamper establishment
};

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

// DATA:

// arrays of parameters:
rate _inflt2attempt_rate[MAX_N_INFLT];
probunit _inflt2delta_probunit[MAX_N_INFLT];
entity_type _e2et[MAX_N_E];

// derived constants:
entity max_e = -1; // largest entity id in use
set<entity> es;
unordered_map<entity_type_pair, set<relationship_or_action_type>> ets2relations;  // possible relations
//unordered_map<event_type, probunit> evt2base_probunit;  // basic success probunit

// variable data:

timepoint current_t = 0;
event current_ev = {};
event_data* current_evd_ = NULL;
// network state:
unordered_map<entity_type, vector<entity>> et2es = {};  // kept to equal inverse of v2vt
unordered_map<entity, leg_set> e2outs = {};
unordered_map<entity, leg_set> e2ins = {};
int n_links = 0; // total no. of current (non-id.) links incl. inverse relationships
// event data:
unordered_map<event, event_data> ev2data = {};
map<timepoint, event> t2be = {};  // kept to equal inverse of ev2data.t

ofstream gexf;
unordered_map<link, timepoint> gexf_edge2start = {};

// I/O:

ostream& operator<< (ostream& os, const link_type& lt) {
    os << et2label[lt.et1]
        << " " << rat2label[lt.rat13] << " "
        << et2label[lt.et3];
    return os;
}
ostream& operator<< (ostream& os, const event_type& evt) {
    os << ec2label[evt.ec]
            << " \"" << et2label[evt.et1]
            << " " << rat2label[evt.rat13] << " "
            << et2label[evt.et3]
            << "\"";
    return os;
}
ostream& operator<< (ostream& os, const event& ev) {
    if (ev.e1 < 0) {
        os << ec2label[ev.ec]
                << " \"some " << et2label[-ev.e1]
                << " " << rat2label[ev.rat13] << " some "
                << et2label[-ev.e3]
                << "\"";
    } else {
        os << ec2label[ev.ec]
                << " \"" << e2label[ev.e1]
                << " " << rat2label[ev.rat13] << " "
                << e2label[ev.e3]
                << "\"";
    }
    return os;
}

/*
void read_config () {
    YAML::Node config = YAML::LoadFile("tricl.config.yaml");

//    for (auto& etspec : config["entity types"]){
//      std::cout << "entity type:" << etspec << endl;
//    }

    // const std::string username = config["username"].as<std::string>();
}
*/

void init_gexf () {
    gexf.open(gexf_filename);
    gexf << R"V0G0N(<?xml version="1.0" encoding="UTF-8"?>
<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">
    <meta>
        <creator>tricl</creator>
        <description>dynamic graph generated by tricl model</description>
    </meta>
    <graph mode="dynamic" defaultedgetype="directed">
        <attributes class="node">
            <attribute id="0" title="entity type" type="string"/>
        </attributes>
        <attributes class="edge">
            <attribute id="1" title="relationship_or_action_type" type="string"/>
        </attributes>
        <nodes>
)V0G0N";
    for (auto& e : es) {
        gexf << "<node id=\"" << e << "\" label=\"" << e2label[e]
             << "\" start=\"0.0\" end=\"" << max_t
             << "\"><attvalues><attvalue for=\"0\" value=\""
             << et2label[_e2et[E(e)]] << "\"/></attvalues></node>" ;
    }
    gexf << R"V0G0N(
        </nodes>
        <edges>
)V0G0N";
}

void dump_data () // dump ev2data to console for debugging
{
    cout << "!" << endl;
    for (auto& [ev2, evd2] : ev2data)
        cout << ev2 << " " << evd2.attempt_rate << " " << evd2.t << endl;
}

void gexf_output_edge (link& l) {
    entity e1 = l.e1, e3 = l.e3; relationship_or_action_type rat13 = l.rat13;
    if (rat13 != RT_ID) {
        gexf << "\t\t\t<edge id=\"" << e1 << "_" << rat13 << "_" << e3 << "_" << current_t
             << "\" source=\"" << e1 << "\" target=\"" << e3
             << "\" start=\"" << gexf_edge2start.at(l) << "\" end=\"" << current_t
             << "\"><attvalues><attvalue for=\"1\" value=\"" << rat2label[rat13]
             << "\"/></attvalues></edge>" << endl;
    }
    gexf_edge2start.erase(l);
}

void finish_gexf () {
    current_t = max_t;
    for (auto& [e1, outs] : e2outs) {
        for (auto& [e3, rat13] : outs) {
            link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 };
            gexf_output_edge(l);
        }
    }
    gexf << R"V0G0N(
        </edges>
    </graph>
</gexf>
)V0G0N";
    gexf.close();
}

// MATH:

// random generators:
random_device ran_dev;
mt19937 random_variable(ran_dev());
uniform_real_distribution<> uniform(0, 1);
exponential_distribution<> exponential(1);

rate effective_rate (rate attempt_rate, probunit success_probunit)
{
    rate r = attempt_rate * probunit2probability(success_probunit);
    assert (r >= 0);
    return r;
}

// LOGICS:

void init_data ()
{
    // init with zeroes:
    for (int i=0; i<MAX_N_INFLT; i++) {
        _inflt2delta_probunit[i] = _inflt2attempt_rate[i] = 0;
    }
    // store actual values:
    for (auto& [inflt, ar] : inflt2attempt_rate) {
        assert (!( inflt.evt.ec == EC_EST && ( inflt.at.rat12 == NO_RAT || inflt.at.rat23 == NO_RAT ) && !(inflt.at == NO_ANGLE)));
        _inflt2attempt_rate[INFLT(inflt)] = ar;
    }
    for (auto& [inflt, spu] : inflt2delta_probunit) {
        assert (!( inflt.evt.ec == EC_EST && ( inflt.at.rat12 == NO_RAT || inflt.at.rat23 == NO_RAT ) ));
        _inflt2delta_probunit[INFLT(inflt)] = spu;
    }
    for (auto& [e, et] : e2et) {
        _e2et[E(e)] = et;
    }
}
void add_entity (entity_type et)
{
    auto e = ++max_e;
    es.insert(e);
    _e2et[E(e)] = et;
    et2es[et].push_back(e);
    e2label[e] = to_string(e);
    e2outs[e] = e2ins[e] = { { .e = e, .r = RT_ID } };
}
entity random_entity (entity_type et)
{
    int pos = floor(uniform(random_variable) * et2es[et].size());
    auto e = et2es[et][pos];
    return e;
}

bool link_exists (link& l)
{
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;
    return (e2outs[e1].count({ .e = e3, .r = rat13}) > 0);
}
void add_link (link& l)
{
    assert ((!link_exists(l)) && ((l.e1 != l.e3) || (l.rat13 == RT_ID)));
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;
    e2outs[e1].insert({ .e = e3, .r = rat13 });
    e2ins[e3].insert({ .e = e1, .r = rat13 });
    if (rat13 != RT_ID) gexf_edge2start[l] = current_t;
    n_links++;
}
void del_link (link& l)
{
    assert (link_exists(l));
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;
    e2outs[e1].erase({ .e = e3, .r = rat13 });
    e2ins[e3].erase({ .e = e1, .r = rat13 });
    if (rat13 != RT_ID) gexf_output_edge(l);
    n_links--;
}

// adapted from set_intersection template:
angles leg_intersection (leg_set& out, leg_set& in)
{
  auto out_it = out.begin(), in_it = in.begin();
  angles result(out.size() + in.size()); // allocate max. size of result
  auto result_it = result.begin();
  while ((out_it != out.end()) && (in_it != in.end()))
  {
      if (out_it->e < in_it->e) {
          ++out_it;
      } else if (in_it->e < out_it->e) {
          ++in_it;
      } else {
          *result_it = { .rat12 = out_it->r, .e2 = out_it->e, .rat23 = in_it->r };
          if (debug) cout << "      leg intersection: " << rat2label[out_it->r] << " " << e2label[out_it->e] << " " << rat2label[in_it->r] << endl;
          ++result_it; ++out_it; ++in_it;
      }
  }
  result.resize(result_it - result.begin()); // shorten result to actual length
  return result;
}

bool event_is_scheduled (event& ev, event_data* evd_)
{
    assert(evd_ == &ev2data[ev]);
    return (evd_->t > -INFINITY);
}
void _schedule_event (event& ev, event_data* evd_)
{
    assert(evd_ == &ev2data[ev]);
    timepoint t;
    auto ar = evd_->attempt_rate;
    assert (ar > 0);
    if (ev.e1 >= 0) { // particular event: use effective rate:
        auto spu = evd_->success_probunit;
        if (ar < INFINITY) {
            t = current_t + exponential(random_variable) / effective_rate(ar, spu);
        } else {
            if (spu > -INFINITY) {
                // event should happen "right away". to make sure all those events
                // occur in random order, we formally schedule them at some
                // random "past" time istead:
                t = current_t - uniform(random_variable);
            } else {
                t = INFINITY;
            }
        }
    } else { // summary event: use only attempt rate:
        t = current_t + exponential(random_variable) / ar;
    }
    t2be[t] = ev;
    evd_->t = t;
    if (verbose) cout << "      scheduled " << ev << " for t=" << t << endl;
}
void schedule_event (event& ev, event_data* evd_)
{
    assert(evd_ == &ev2data[ev]);
    assert (!event_is_scheduled(ev, evd_));
    _schedule_event(ev, evd_);
}
void reschedule_event (event& ev, event_data* evd_)
{
    assert(evd_ == &ev2data[ev]);
    assert (event_is_scheduled(ev, evd_));
    t2be.erase(evd_->t);
    _schedule_event(ev, evd_);
}
void add_event (event& ev)
{
    auto ec = ev.ec; auto e1 = ev.e1, e3 = ev.e3; auto rat13 = ev.rat13;
    assert ((rat13 != RT_ID) && (e1 != e3));
    auto et1 = _e2et[E(e1)], et3 = _e2et[E(e3)];
    event_type evt = { .ec = ec, .et1 = et1, .rat13 = rat13, .et3 = et3 };
    if (evt2base_probunit.count(evt) > 0) { // event can happen at all:
        if (debug) cout << "     adding event: " << ev << endl;
        // find and store attempt rate and success probunit by looping through adjacent legs and angles:
        rate ar = 0;
        auto spu = evt2base_probunit[evt];
        auto outs1 = e2outs[e1], ins3 = e2ins[e3];
        // legs:
        if (ec == EC_TERM) {
            for (auto& [e2, rat12] : outs1) {
                influence_type inflt = { .evt = evt, .at = { .rat12 = rat12, .et2 = _e2et[E(e2)], .rat23 = NO_RAT } };
                ar += _inflt2attempt_rate[INFLT(inflt)];
                spu += _inflt2delta_probunit[INFLT(inflt)];
            }
            for (auto& [e2, rat23] : ins3) {
                influence_type inflt = { .evt = evt, .at = { .rat12 = NO_RAT, .et2 = _e2et[E(e2)], .rat23 = rat23 } };
                ar += _inflt2attempt_rate[INFLT(inflt)];
                spu += _inflt2delta_probunit[INFLT(inflt)];
            }
        }
        // angles:
        int na = 0; // number of influencing angles
        angles as = leg_intersection(outs1, ins3);
        for (auto a_it = as.begin(); a_it < as.end(); a_it++) {
            influence_type inflt = { .evt = evt, .at = { .rat12 = a_it->rat12, .et2 = _e2et[E(a_it->e2)], .rat23 = a_it->rat23 } };
            if (debug) cout << "      influences of angle " << rat2label[a_it->rat12] << " " << e2label[a_it->e2] << " " << rat2label[a_it->rat23] << ":" << endl;
            auto dar = _inflt2attempt_rate[INFLT(inflt)];
            auto dsl = _inflt2delta_probunit[INFLT(inflt)];
            if (dar != 0 || dsl != 0) { // angle can influence event
                na++;
                if (debug) {
                    if (dar != 0) cout << "       on attempt rate:" << dar << endl;
                    if (dsl != 0) cout << "       on success probunit:" << dsl << endl;
                }
                ar += dar;
                spu += dsl;
            } else if (debug) cout << "       none" << endl;
        }
        // only add a non-termination event individually if at least one angle existed
        // (pure spontaneous non-termination events are handled summarily via entity types to keep maps sparse):
        if ((ec == EC_TERM) || (na > 0)) {
            // add spontaneous event probs.:
            influence_type inflt = { .evt = evt, .at = NO_ANGLE };
            ar += _inflt2attempt_rate[INFLT(inflt)];
            spu += _inflt2delta_probunit[INFLT(inflt)];
            ev2data[ev] = { .n_angles = na, .attempt_rate = ar, .success_probunit = spu, .t = -INFINITY };
            if (debug) cout << "      attempt rate " << ar << ", success prob. " << probunit2probability(spu) << endl;
            schedule_event(ev, &ev2data[ev]);
        }
    } else {
        if (debug) cout << "     not adding impossible event: " << ev << endl;
    }
}
void remove_event (event& ev, event_data* evd_)
{
    assert (event_is_scheduled(ev, evd_));
    ev2data.erase(ev);
    t2be.erase(evd_->t);
}

void add_or_delete_angle (event_class ec, // add or delete?
        entity e1, entity_type et1, relationship_or_action_type rat12,
        entity e2, entity_type et2, relationship_or_action_type rat23,
        entity e3, entity_type et3)
{
    if (debug) cout << "   " << ec2label[ec] << " angle \"" << e2label[e1] << " " << rat2label[rat12] << " "
            << e2label[e2] << " " << rat2label[rat23] << " " << e2label[e3] << "\"" << endl;
    for (auto& rat13 : ets2relations[{ et1, et3 }]) {
        event_class ec13 = (e2outs[e1].count({ .e = e3, .r = rat13 }) == 0) ? EC_EST : EC_TERM;
        event ev = { .ec = ec13, .e1 = e1, .rat13 = rat13, .e3 = e3 };
        auto evd_ = &ev2data[ev];

        // FIXME: deal with case where event is not scheduled yet!!!

        event_type evt = { .ec = ec13, .et1 = et1, .rat13 = rat13, .et3 = et3 };
        if (debug) cout << "    possibly updating event: " << ec2label[ec13] <<  " \"" << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << "\"" << endl;
        if (evt2base_probunit.count(evt) > 0) {
            if (debug) cout << "     event type " << evt << " has a base success prob. of " << probunit2probability(evt2base_probunit[evt]) << endl;
            influence_type inflt = {
                    .evt = evt,
                    .at = { .rat12 = rat12, .et2 = et2, .rat23 = rat23 }
            };
            auto dar = _inflt2attempt_rate[INFLT(inflt)], dsl = _inflt2delta_probunit[INFLT(inflt)];
            if (dar != 0 || dsl != 0) {
                if (debug) cout << "     angle may influence attempt or success" << endl;
                if (ec == EC_EST) { // angle is added:
                    evd_->n_angles += 1;
                    if (dar != 0) {
                        evd_->attempt_rate += dar;
                    }
                    if (dsl != 0) {
                        evd_->success_probunit += dsl;
                    }
                    auto& told = evd_->t;
                    if (told > -INFINITY) t2be.erase(told);
                    schedule_event(ev, evd_);
                } else { // angle is removed:
                    evd_->n_angles -= 1;
                    if (dar != 0) {
                        evd_->attempt_rate -= dar;
                    }
                    if (dsl != 0) {
                        evd_->success_probunit -= dsl;
                    }
                    if ((ec13 != EC_TERM) && (evd_->n_angles == 0)) { // only spontaneous non-termination event is left:
                        // remove specific event:
                        remove_event(ev, evd_); // FIXME: sometimes event was not scheduled before...
                    } else {
                        auto& told = evd_->t;
                        if (told > -INFINITY) t2be.erase(told);
                        schedule_event(ev, evd_);
                    }
                }
            } else {
                if (debug) cout << "     but angle may not influence attempt or success" << endl;
            }
        }
    }
}

void init_entities ()
{
    // inspect pre-registered entities:
    auto vt2remaining_n = et2n;
    for (auto& [e, et] : e2et) {
        assert(e >= 0);
        es.insert(e);
        if (e > max_e) max_e = e;
        if (e2label.count(e) == 0) e2label[e] = to_string(e);
        if (e2outs.count(e) == 0) e2outs[e] = {};
        if (e2ins.count(e) == 0) e2ins[e] = {};
        if (et2es.count(et) == 0) et2es[et] = {};
        et2es[et].push_back(e);
        vt2remaining_n[et]--;
    }
    // generate remaining entities:
    for (auto& [et, n] : vt2remaining_n) {
        cout << "entity type \"" << et2label[et] << "\" has " << et2n[et] << " entities" << endl;
        if (n < et2n[et]) {
            cout << " of which were preregistered:" << endl;
            for (auto& e : et2es[et]) cout << "  " << e2label[e] << endl;
        }
        assert (n >= 0);
        while (n > 0) {
            add_entity(et);
            n--;
        }
        assert ((entity) et2es[et].size() == et2n[et]);
    }
}

void init_relationship_or_action_types ()
{
    // verify symmetry of relationship inversion map:
    cout << "relationship types (with inverses) and action types:" << endl;
    assert (rat2inv[RT_ID] == RT_ID);
    for (auto& [r, la] : rat2label) {
        if (r_is_action_type[r]) {
            cout << " action type \"" << la << "\"" << endl;
        } else {
            auto inv = rat2inv[r];
            if (inv != NO_RAT) {
                assert (rat2inv[inv] == r);
                cout << " relationship type \"" << la << "\" (inverse: \"" << rat2label[inv] << "\")" << endl;
            } else {
                cout << " relationship type \"" << la << "\"" << endl;
            }
        }
    }
    // register possible relationship types by entity type pair, and compute probunits:
    cout << "base success probability units:" << endl;
    for (auto& [evt, pu] : evt2base_probunit) {
        auto rat13 = evt.rat13;
        assert (rat13 != RT_ID);
        entity_type_pair ets = { evt.et1, evt.et3 };
        if (ets2relations.count(ets) == 0) ets2relations[ets] = {};
        ets2relations[ets].insert(rat13);
        cout << " " << evt << ": " << pu << endl;
    }
}

void init_generic_events ()
{
    cout << "initial scheduling of generic events" << endl;
    // summary events for purely spontaneous additions:
    for (auto& [ets, relations] : ets2relations) {
        auto et1 = ets.et1, et3 = ets.et3;
        for (auto& rat13 : relations) {
            event ev = { .ec = EC_EST,
                    .e1 = (entity)-et1, // in spontaneous events, fields e1 and e3 are used to store entity types with negative sign
                    .rat13 = rat13,
                    .e3 = (entity)-et3 };
            event_type evt = { .ec = EC_EST, .et1 = et1, .rat13 = rat13, .et3 = et3 };
            influence_type inflt = { .evt = evt, .at = NO_ANGLE };
            auto ar = _inflt2attempt_rate[INFLT(inflt)];
            if (ar > 0) {
                ev2data[ev] = { .n_angles = 0, .attempt_rate = ar * et2n[et1] * et2n[et3], .success_probunit = evt2base_probunit[evt] + _inflt2delta_probunit[INFLT(inflt)]};
                schedule_event(ev, &ev2data[ev]);
            }
        }
    }
}

bool pop_next_event ()
{
    bool found = false;
    // find next event:
    while ((!found) && (current_t < max_t)) {
        auto tev = t2be.begin();
        if (tev == t2be.end()) {
            current_t = max_t;
            return false;
        }
        timepoint t = tev->first;
        if (t >= max_t) {
            current_t = max_t;
            return false;
        }
        event ev = tev->second;
        if (t > current_t) {
            current_t = t;
        } // otherwise it's an event happening "right now" scheduled formally for a past time to ensure a random order of those events.
        if (ev.e1 < 0) { // event is purely spontaneous addition attempt, so has only types specified
            event summary_ev = ev;
            // draw actual entities at random from given types:
            entity_type et1 = -summary_ev.e1, et3 = -summary_ev.e3;
            auto e1 = random_entity(et1), e3 = random_entity(et3);
            auto rat13 = summary_ev.rat13;
            link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 };
            if ((e1 != e3) && (!link_exists(l))) {
                // link does not already exist
                event actual_ev = { .ec = EC_EST, .e1 = e1, .rat13 = rat13, .e3 = e3 };
                if (ev2data.count(actual_ev) == 0) {
                    // no non-pure addition event for this pair was timed,
                    // so check if attempt succeeds:
                    event_type evt = { .ec = EC_EST, .et1 = et1, .rat13 = rat13, .et3 = et3 };
                    auto spu = evt2base_probunit.at(evt);
                    for (auto& [e2, rat12] : e2outs[e1]) {
                        influence_type inflt = { .evt = evt, .at = { .rat12 = rat12, .et2 = _e2et[E(e2)], .rat23 = NO_RAT } };
                        if (inflt2delta_probunit.count(inflt) > 0) spu += inflt2delta_probunit.at(inflt);
                    }
                    for (auto& [e2, rat23] : e2ins[e3]) {
                        influence_type inflt = { .evt = evt, .at = { .rat12 = NO_RAT, .et2 = _e2et[E(e2)], .rat23 = rat23 } };
                        if (inflt2delta_probunit.count(inflt) > 0) spu += inflt2delta_probunit.at(inflt);
                    }
                    if (uniform(random_variable) < probunit2probability(spu)) { // success
                        // update simulation time and return event:
                        current_ev = actual_ev;
                        found = true;
                    } else {
                        if (verbose) cout << "at t=" << current_t << ", link addition failed" << endl;
                    }
                }
            } else { // otherwise try again...
                if (verbose) cout << "at t=" << current_t << ", link to add existed already" << endl;
            }
            // set next_occurrence of this summary event:
            reschedule_event(summary_ev, &ev2data[summary_ev]);
        } else { // event has specific entities
            // update simulation time and return event:
            current_ev = ev;
            auto evd_ = current_evd_ = &ev2data[ev];
            // remove it from all relevant data:
            remove_event(ev, evd_);
            found = true;
        }
    }
    return found;
}

#define IF_12 if (which == 0)
#define IF_23 if (which == 1)
#define EPS 1e-20;

void update_adjacent_events (event& ev)
{
    auto ec_ab = ev.ec; auto ea = ev.e1, eb = ev.e3; auto rab = ev.rat13;
    entity e1, e2, e3; entity_type et1, et2, et3; relationship_or_action_type rat12, rat23;
    leg_set legs;
    if (verbose) cout << "  updating adjacent events of " << ev << endl;
    // loop through all adjacent events e1->e3:
    for (int which = 0; which < 2; which++) {
        IF_12 { // ea->eb plays role of e1->e2:
            if (verbose) cout << "   angles with this as 1st leg:" << endl;
            e1 = ea; rat12 = rab; e2 = eb;
            et1 = _e2et[E(e1)]; et2 = _e2et[E(e2)];
            legs = e2outs[eb]; // rat23, e3
        }
        IF_23 { // ea->eb plays role of e2->e3:
            if (verbose) cout << "   angles with this as 2nd leg:" << endl;
            legs = e2ins[ea]; // e1, rat12
            e2 = ea; rat23 = rab; e3 = eb;
            et2 = _e2et[E(e2)]; et3 = _e2et[E(e3)];
        }
        for (auto& l : legs) {
            IF_12 {
                rat23 = l.r; e3 = l.e; et3 = _e2et[E(e3)];
            }
            IF_23 {
                e1 = l.e; rat12 = l.r; et1 = _e2et[E(e1)];
            }
            if (e1 != e3) add_or_delete_angle(ec_ab, e1, et1, rat12, e2, et2, rat23, e3, et3);
        }
    }
}

void add_reverse_event (event& old_ev)
{
    auto e1 = old_ev.e1, e3 = old_ev.e3; auto rat13 = old_ev.rat13;
    event_class ec = (old_ev.ec == EC_TERM) ? EC_EST : EC_TERM;
    event ev = { .ec = ec, .e1 = e1, .rat13 = rat13, .e3 = e3 };
    if (debug) cout << "    adding reverse event: " << ev << endl;
    add_event(ev);
}

void perform_event (event& ev)
{
    if (debug) cout << " performing event: " << ev << endl;
    auto ec = ev.ec; auto e1 = ev.e1, e3 = ev.e3; auto rat13 = ev.rat13, r31 = rat2inv[rat13];
    event companion_ev = { .ec = ec, .e1 = e3, .rat13 = r31, .e3 = e1 };
    link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 }, inv_l = { .e1 = e3, .rat13 = r31, .e3 = e1 };
    // add or remove link:
    if (ec == EC_EST) {
        add_link(l);
    } else {
        del_link(l);
    }
    add_reverse_event(ev);
    update_adjacent_events(ev);
    // also perform companion event that affects inverse link:
    if (r31 != NO_RAT) {
        auto companion_evd_ = &ev2data[companion_ev];
        if (companion_evd_->t > -INFINITY) {
            remove_event(companion_ev, companion_evd_);
            if (debug) cout << " unscheduling companion event" << endl;
        }
        if (ec == EC_EST) {
            if (verbose) cout << " performing companion event: adding inverse link \"" << e2label[e3]
                << " " << rat2label[r31] << " " << e2label[e1] << "\"" << endl;
            add_link(inv_l);
        } else {
            if (verbose) cout << " performing companion event: deleting inverse link \"" << e2label[e3]
                << " " << rat2label[r31] << " " << e2label[e1] << "\"" << endl;
            del_link(inv_l);
        }
        add_reverse_event(companion_ev);
        update_adjacent_events(companion_ev);
    }
}

void do_random_link (probability p, entity_type e1, relationship_or_action_type rat13, entity_type e3) {
    if (uniform(random_variable) < p) {
        link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 };
        if (!link_exists(l)) {
            event ev = { .ec = EC_EST, .e1 = e1, .rat13 = rat13, .e3 = e3 };
            auto evd_ = &ev2data[ev];
            if (evd_->t > -INFINITY) remove_event(ev, evd_);
            perform_event(ev);
        }
    }
}

void init_links ()
{
    cout << "performing events that add initial links" << endl;
    // identity relationship:
    for (auto& e : es) {
        leg l = { .e = e, .r = RT_ID };
        e2outs[e].insert(l);
        e2ins[e].insert(l);
    }
    // preregistered links:
    for (auto l : initial_links) {
        assert (l.rat13 != RT_ID);
        event ev = { .ec = EC_EST, .e1 = l.e1, .rat13 = l.rat13, .e3 = l.e3 };
        auto evd_ = &ev2data[ev];
        if (!link_exists(l)) {
            if (evd_->t > -INFINITY) remove_event(ev, evd_);
            perform_event(ev);
        }
    }
    // random links:
    // block model:
    unordered_map<entity, int> e2block = { };
    for (auto& e : es) {
        e2block[e] = floor(uniform(random_variable) * et2n_blocks[_e2et[E(e)]]);
    }
    for (auto& [lt, pw] : lt2initial_prob_within) { // TODO: add lt2initial_prob_between
        if (pw > 0) {
            assert (lt.rat13 != RT_ID);
            auto et1 = lt.et1, et3 = lt.et3; auto rat13 = lt.rat13;
            for (auto& e1 : et2es[et1]) {
                for (auto& e3 : et2es[et3]) {
                    if ((e3 != e1) && ((e3 > e1) || (rat2inv[rat13] != rat13))) {
                        float p = (e2block[e1] == e2block[e3]) ? lt2initial_prob_within[lt] : lt2initial_prob_between[lt];
                        do_random_link(p, e1, rat13, e3);
                    }
                }
            }
        }
    }
    // random geometric model:
    unordered_map<entity, vector<float>> e2coords = { };
    for (auto& [et, dim] : et2dim) {
        for (auto& e : et2es[et]) {
            for (int d=0; d<dim; d++) {
                e2coords[e].push_back(uniform(random_variable));
            }
        }
    }
    for (auto& [lt, ex] : lt2spatial_decay) {
        assert (lt.rat13 != RT_ID);
        auto et1 = lt.et1, et3 = lt.et3; auto rat13 = lt.rat13; auto dim = et2dim[et1];
        assert (et2dim[et3] == dim);
        for (auto& e1 : et2es[et1]) {
            for (auto& e3 : et2es[et3]) {
                if ((e3 != e1) && ((e3 > e1) || (rat2inv[rat13] != rat13))) {
                    float dist2 = 0.0;
                    for (int d=0; d<dim; d++) {
                        dist2 += pow(e2coords[e1][d] - e2coords[e3][d], 2);
                    }
                    float p = exp(-ex * sqrt(dist2));
                    do_random_link(p, e1, rat13, e3);
                }
            }
        }
    }
}

void init ()
{
    cout << "INITIALIZING..." << endl;
    init_data();
    init_entities();
    init_relationship_or_action_types();
    init_generic_events();
    init_links();
    init_gexf();
    cout << "...INITIALIZATION FINISHED." << endl << endl;
}

void step ()
{
    if (pop_next_event()) {
        if (!quiet) {
            cout << "at t=" << current_t << " " << current_ev << ", density " << ((double)n_links)/(max_e*max_e) << endl;
        }
        perform_event(current_ev);
        if (quiet) {
            cout << "                          \r" << round(current_t) << ": " << ((double)n_links)/(max_e*max_e) << "  ";
        }
        if (verbose) cout << " " << t2be.size() << " events on stack" << endl;
    }
    else {
        current_t = max_t;
    }
}

void finish ()
{
    finish_gexf();
}

int main ()
{
//    read_config();
    init();
    while (current_t < max_t) step();
    if (verbose) {
        cout << "\nat t=" << current_t << ", " << t2be.size() << " events on stack: " << endl;
        for (auto& [t, ev] : t2be) {
            cout << " " << ev << " at " << t << endl;
        }
    }
    finish();
    cout << endl;
    return 0;
}

