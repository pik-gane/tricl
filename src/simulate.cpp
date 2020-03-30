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
