/** TriCl, a generic model of social dynamics
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Oct 18, 2019
 *
 * @file
 */

// DEBUG or RELEASE mode?
// for debug mode comment the following line, for release mode leave it uncommented:
//#define NDEBUG

/** TODO:
 *
 * model capabilities:
 * - allow legs to influence success prob. of establishment
 * - allow legs to attempt establishment
 * - add acts
 *
 * simulation options:
 * - add command line setting of metaparameters
 * - as an alternative to tmax and max_n_events, add max_wall
 * - add particle filtering mode
 *
 * network theory stuff:
 * - output triangle counts clustering coefficient(s)
 * - add more random network models
 *
 * convenience:
 * - improve quiet status output
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
 * event: an event class plus a link, e.g.
 *        - "establishment of the link 'John - is a subscriber of - the BBC'" (event class EC_EST)
 *        - "termination of the link 'John - is friends with - Alice'" (event class EC_TERM)
 *        - "occurrence of the act 'Alice - utters that - Elvis lives'" (event class EC_ACT)
 *
 * leg: a relationship or action type rat plus an entity.
 *      a leg can influence an adjacent termination event but no establishment or act occurrence events
 *
 * angle: a middle entity plus pair of relationship or action types.
 *        an angle can influence any adjacent event.
 *        an angle can also encode a leg if either relationship or action type equals the special value NO_RT
 *
 * influence: a event plus a leg or angle that influences it (or the special value NO_ANGLE if the event happens spontaneously)
 *
 * attempt rate of an influence: the probability rate with which the influence attempts the corresponding event
 *
 * success probability units of an event: probability that an attempted event actually happens, transformed via a sigmoidal function.
 *
 * the following diagrams show the relationships between these types of things:
 *
 *       e1 –––––––––rat13–––––––––> e3
 *     entity     relationship     entity
 *       .       or action type       .
 *       .             .              .
 *       .       \______________________/
 *       .           an out-leg for e1
 *       .             .              .
 *     \______________________/       .
 *         an in-leg for e3           .
 *       .             .              .
 *     \________________________________/
 *                   a link
 *
 *
 *       e1 ––rat12––> e2 ––rat23––> e3
 *      \______________________________/
 *                  an angle
 */

/** ABBREVIATIONS used in variable naming.
 *
 * e: entity
 * e1: source entity of a link or event
 * e2: middle entity of an angle or influence, or "other" entity of a leg
 * e3: target entity of a link or event
 *
 * et: entity type
 * et1: source entity type of a link or event type
 * et2: middle entity type of an angle or influence type, or other entity type of a leg type
 * et3: target entity type of a link or event type
 *
 * ev: event
 * evt: event type
 *
 * infl: influence
 * inflt: influence type
 *
 * l: link
 * lt: link type
 *
 * rat: relationship or action type
 * rat13: type of rel. or action from source to target of a link or event
 * rat12: type of rel. or action from source to middle/other entity
 * rat23: type of rel. or action from middle/other entity to target
 *
 * t: timepoint
 *
 * trailing _: pointer
 * trailing _it: pointer used as iterator
 *
 * x2y: map mapping xs to ys (typically a map, unordered_map, or vector)
 *
 * OTHER NAMING CONVENTIONS:
 *
 * void functions are named by imperatives (e.g., "do_this")
 * non-void functions are named by nouns (e.g., "most_important_thing")
 *
 */

// INCLUDES:

#include "global_variables.h"
#include "io.h"
#include "config.h"
#include "init.h"
#include "simulate.h"
#include "finish.h"
#include "debugging.h"
//#include "pfilter.h"

string config_yaml_filename; ///< filename of configuration file

/** main function of the tricl executable */
int main (int argc, char *argv[]) {
    try
    { // because our code sometimes deliberately throws strings upon errors

        // stuff before simulation:
        read_config(argc, argv);
        init();
        if (debug) verify_data_consistency();

        // actual simulation:
        while (true)
        {
            if (!step()) break;
        }

        // stuff after simulation:
        finish();
    }
    catch (const char* msg)
    {
      cerr << "ERROR: exiting with message: " << msg << endl;
    }
    return 0;
}
