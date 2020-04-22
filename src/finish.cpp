/** Do stuff at end of the simulation.
 *
 *  \file
 */

#include "global_variables.h"
#include "debugging.h"
#include "io.h"
#include "gexf.h"
#include "finish.h"

/** Do stuff at end of the simulation.
 */
void finish ()
{
    // forward to end of simulation time:
    current_t = max_t;  // TODO: do we need this?

    if (verbose) {
        cout << "\nat t=" << current_t << ", " << t2ev.size() << " events on stack: " << endl;
        for (auto& [t, ev] : t2ev) {
            cout << " " << ev << " at " << t << endl;
        }
    }
    log_state();
    cout << endl;

    finish_gexf();

    if (debug) verify_data_consistency();
}
