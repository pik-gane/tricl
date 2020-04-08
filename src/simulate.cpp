/*
 * simulate.cpp
 *
 *  Created on: Oct 18, 2019
 *      Author: Jobst Heitzig, PIK
 */

#include <iostream>

#include "global_variables.h"
#include "event.h"
#include "simulate.h"

using namespace std;

bool step ()
{
    if ((n_events < max_n_events) && pop_next_event()) {
        ++n_events;
        perform_event(current_ev);
        if (quiet) {
            cout << round(current_t) << ": " << ((double)n_links)/((double)max_e*(double)max_e) << "                            \r";
        }
        if (debug) cout << " " << t2be.size() << " events on stack" << endl << endl;
        return true;
    } else {
        return false;
    }
}
