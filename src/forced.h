// make sure this file is only included once:
#ifndef INC_FORCED_H
#define INC_FORCED_H

/** Forced (exogenous) events and other events at deterministic times.
 *
 *  \file
 *
 *  Forced events are read from a csv file (option --forced-in, same columns as event files) and applied at their
 *  times in simulation and replay. They must be events of link types without dynamics in the config file
 *  (e.g. a timetable encoded as links of all agents to hub entities), are not written to the events output, and
 *  add nothing to the log-likelihood, which is therefore conditional on them.
 *  Together with the expiries of actions (see actions.h), they are the events at deterministic times, which are
 *  processed in time order whenever model time advances.
 */

#include "global_variables.h"

extern string forced_in_filename;   ///< Name of (or path to) the csv file of forced events (if "", none)
extern bool performing_forced;      ///< Whether the event currently performed is a forced event

void read_forced_events ();
bool any_deterministic_events ();
timepoint next_deterministic_time ();
void process_deterministic_until (timepoint t);

#endif
