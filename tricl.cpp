/*
 * tricl.cpp
 *
 *  Created on: Oct 18, 2019
 *      Author: Jobst Heitzig, PIK
 */

/* USAGE:
 *
 * CONFIGURE by editing parameter_values.cpp
 *
 * COMPILE with:
 *   cd build/default
 *   cmake ../../
 *   cmake --build .
 *
 * RUN with ./tricl
 *
 * PROFILE with valgrind --tool=callgrind --callgrind-out-file=/tmp/callgrind.out ./tricl; kcachegrind /tmp/callgrind.out
 */

/* TODO:
 * verify CMake build process according to https://cmake.org/cmake/help/latest/guide/tutorial/index.html
 * clarify how to deal with "assert" in production mode (all source files!)
 * add event preclasses EC_A_EST, EC_A_ACT for attempts of generic events, then allow legs to influence their success prob.
 * use const args for optimization
 * take into account influence of legs on addition!
 * add acts
 * read parameters from YAML file and add --help message
 * output clustering coefficient(s)
 * add particle filtering mode
 * replace inflt2attempt_rate by evt2base_attempt_rate and inflt2delta_attempt_rate
 * add control over random seed
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

#include <iostream>

// local includes:

#include "data_model.h"
#include "global_variables.h"
#include "probability.h"

#include "entity.h"
#include "link.h"
#include "angle.h"
#include "event.h"

//#include "config.h"
#include "io.h"
#include "gexf.h"
#include "init.h"
#include "finish.h"

//#include "pfilter.h"

using namespace std;

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

