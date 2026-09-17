/** Run the actual simulation.
 *
 *  \file
 */

#include "global_variables.h"
#include "event.h"
#include "debugging.h"
#include "simulate.h"

#define RECOMPUTE_TOTAL_ER_EVERY 65536  ///< No. of events after which the total effective rate is recomputed exactly to remove floating point drift

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
        if (n_events % RECOMPUTE_TOTAL_ER_EVERY == 0) {
            // remove floating point drift from the incrementally maintained total effective rate:
            rate exact = compute_total_finite_er();
            if (fabs(exact - total_finite_effective_rate) > 1e-6 * max(1.0, fabs(exact)))
                cerr << "WARNING: total effective rate " << total_finite_effective_rate << " deviates from exact value " << exact
                     << " after " << n_events << " events. Please report this as a bug." << endl;
            total_finite_effective_rate = max(0.0, exact);
        }
        return true;
    } else {
        return false;
    }
}
