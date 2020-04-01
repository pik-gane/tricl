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

void step ()
{
    if (pop_next_event()) {
        perform_event(current_ev);
        if (quiet) {
            cout << round(current_t) << ": " << ((double)n_links)/(max_e*max_e) << "                            \r";
        }
        if (debug) cout << " " << t2be.size() << " events on stack" << endl << endl;
    }
    else {
        current_t = max_t;
    }
}
