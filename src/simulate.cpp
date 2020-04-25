/** Run the actual simulation.
 *
 *  \file
 */

#include "global_variables.h"
#include "event.h"
#include "simulate.h"

/** Perform next step.
 *
 *  \returns whether there was another step to perform.
 */
bool step ()
{
    if ((n_events < max_n_events) && pop_next_event()) {
        ++n_events;
        perform_event(current_ev, current_evd_);
        if (debug) cout << " " << t2ev.size() << " events on stack" << endl << endl;
        return true;
    } else {
        return false;
    }
}
