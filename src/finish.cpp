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
    // forward to end of simulation time
    // (unless the simulation ended due to the event limit, in which case current_t is the time of the last event,
    // or the time limit is infinite):
    if ((n_events >= max_n_events) || (max_t == INFINITY) || wall_time_exceeded) {
        if (!quiet) cout << (wall_time_exceeded ? "wall-clock time limit reached, " : "") << "stopped at t=" << current_t << " after " << n_events << " events." << endl;
    } else {
        current_t = max_t;
    }

    if (verbose) {
        cout << "\nat t=" << current_t << ", " << ev2data.size() << " events scheduled (total rate " << schedule.total()
             << ", " << immediate_events.size() << " immediate): " << endl;
        for (auto& [ev, evd] : ev2data) {
            cout << " " << ev << " at rate " << evd.effective_rate << endl;
        }
    }
    log_state(true);
    if (!silent) cout << endl;

    write_stats_until(current_t);
    close_stats_out();
    finish_gexf();
    close_links_out();
    close_events_out();

    if (debug) verify_data_consistency();

    if (only_output_logl) cout << cumulative_logl << endl;
    if (output_summary) output_json_summary();
}
