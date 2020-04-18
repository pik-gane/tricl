/*
 * finish.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include "global_variables.h"
#include "debugging.h"
#include "io.h"
#include "gexf.h"
#include "finish.h"

void finish ()
{
    current_t = max_t;
    if (verbose) {
        cout << "\nat t=" << current_t << ", " << t2be.size() << " events on stack: " << endl;
        for (auto& [t, ev] : t2be) {
            cout << " " << ev << " at " << t << endl;
        }
    }
    log_status();
    finish_gexf();
    if (debug) verify_data_consistency();
}
