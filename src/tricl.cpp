/** TriCl, a generic model of social dynamics
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Oct 18, 2019
 *
 * @file
 */

/** USAGE.
 *
 * BUILD:
 *   cd tricl
 *   mkdir -p build/default
 *   cd build/default
 *   cmake ../../
 *   cmake --build .
 *
 * CONFIGURE: copy parameters_TEMPLATE.yaml to some file my_parameters.yaml and edit it
 *
 * RUN: ./tricl my_parameters.yaml
 *
 * PROFILE: valgrind --tool=callgrind --callgrind-out-file=/tmp/callgrind.out ./tricl my_parameters.yaml; kcachegrind /tmp/callgrind.out
 *
 */

// DEBUG or RELEASE mode?
// for debug mode comment the following line, for release mode leave it uncommented it:
//#define NDEBUG

/** TODO:
 *
 * model capabilities:
 * - allow legs to influence success prob. of establishment (add event preclasses EC_A_EST, EC_A_ACT for attempts)
 * - allow legs to attempt establishment
 * - add acts
 *
 * simulation options:
 * - as an alternative to tmax, add max_n_events and max_wall
 * - add control over random seed
 * - allow to read initial state from csv files (one for each link type) or gexf file
 * - add particle filtering mode
 *
 * network theory stuff:
 * - output clustering coefficient(s)
 * - add more random network models
 *
 * convenience:
 * - add --help message and potential parameter overrides on the command line
 * - remove redundant output from init that is already printed in config
 *
 * optimization:
 * - use const args as much as possible in inner loops
 *
 */

/** TERMINOLOGY (see also overleaf document).
 *
 * entity type: a type of concrete or abstract entity (see below),
 *              e.g. "an individual", "a news channel", "a social group", "an opinion", or "an infection state"
 *
 * entity: any concrete or abstract object that can stand in a relationship with other entities,
 *         e.g. "John", "the BBC", "catholics", "Elvis lives", or "infected with Dengue"
 *
 * relationship type: any concrete or abstract type of directed relationship two entities can stand in,
 *                    e.g. "is friends with", "is a subscriber of", "belongs to", "holds", or "is".
 *                    the special relationship type "=" (with id RT_ID) encodes the identity of an entity with itself.
 *
 * action type: any type of thing that can happen at a singular time point between two entities,
 *              e.g. "kisses" or "utters that"
 *
 * act: a pair of entities plus an action type plus a time-point.
 *      we say the act "occurs" at that time.
 *
 * link type: a pair of entity types plus a relationship or action type,
 *            e.g. "an individual - is a subscriber of - a news channel" or "an individual - utters - an opinion"
 *
 * link: a pair of entities plus a relationship or action type,
 *       e.g. "John - is a subscriber of - the BBC".
 *       if it has a relationship type, a link "exists" as long as the entities stand in the respective relationship.
 *       if it has an action type, a link "flashes" whenever a corresponding act occurs.
 *       the source entity id is called "e1" in code, the destination entity id "e3" (because of angles, see below)
 *
 * relationship: a link of relationship type
 *
 * action: a link of action type (not to be confused with an act!)
 *
 * impact of an action: a nonnegative real number that is increased by one whenever the link flashes and decays exponentially at a certain rate.
 *                      a link's impact determines the link's influence on adjacent events.
 *
 * basic event: an event class plus a link, e.g.
 *              - "establishment of the link 'John - is a subscriber of - the BBC'" (event class EC_EST)
 *              - "termination of the link 'John - is friends with - Alice'" (event class EC_TERM)
 *              - "occurrence of the act 'Alice - utters that - Elvis lives'" (event class EC_ACT)
 *
 * leg: a relationship or action type rat plus an entity.
 *      a leg can influence an adjacent termination event but no establishment or act occurrence events
 *
 * angle: a middle entity plus pair of relationship or action types.
 *        an angle can influence any adjacent event.
 *        an angle can also encode a leg if either relationship or action type equals the special value NO_RT
 *
 * influence: a basic event plus a leg or angle that influences it (or the special value NO_ANGLE if the basic event happens spontaneously)
 *
 * attempt rate of an influence: the probability rate with which the influence attempts the corresponding event
 *
 * success probability units of a basic event: probability that an attempted event actually happens, transformed via a sigmoidal function.
 *
 */

/** ABBREVIATIONS used in variable naming.
 *
 * x2y: map mapping xs to ys
 * e[t]: entity (aka node) [type] id
 * ev[t]: event [type] id
 * infl[t]: influence [type] id
 * l[t]: link [type] id
 * rat: relationship or action type id
 * t: timepoint
 * trailing _: pointer
 * trailing _it: pointer used as iterator
 *
 * OTHER NAMING CONVENTIONS:
 *
 * void functions are named by imperatives (e.g., "do_this")
 * non-void functions are named by nouns (e.g., "most_important_thing")
 *
 */

// INCLUDES:

// libraries:

#include <iostream>

// local includes:

#include "debugging.h"
#include "global_variables.h"
#include "config.h"
#include "init.h"
#include "simulate.h"
#include "finish.h"
//#include "pfilter.h"

using namespace std;

string config_yaml_filename; ///< filename of configuration file

/** parse the command line arguments */
void parse_args (int argc, const char *argv[]) {
    if (argc < 2) {
        throw "USAGE: tricl my_config.yaml";
    }
    config_yaml_filename = argv[1];
}

/** main function of the tricl executable */
int main (int argc, const char *argv[]) {
    try {
        parse_args(argc, argv);
        read_config();
        init();
        if (debug) verify_data_consistency();

        while (true) {
            if (!step()) break;
        }

        if (verbose) {
            cout << "\nat t=" << current_t << ", " << t2be.size() << " events on stack: " << endl;
            for (auto& [t, ev] : t2be) {
                cout << " " << ev << " at " << t << endl;
            }
        }
        finish();
        cout << endl;
    } catch (const char* msg) {
      cerr << "ERROR: exiting with message: " << msg << endl;
    }
    return 0;
}
