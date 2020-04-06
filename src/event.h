/*
 * event.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_EVENT_H
#define INC_EVENT_H

#include "data_model.h"

using namespace std;

bool event_is_scheduled (event& ev, event_data* evd_);

void schedule_event (event& ev, event_data* evd_, double left_tail, double right_tail);

void reschedule_event (event& ev, event_data* evd_, double left_tail, double right_tail);

void add_event (event& ev);

void remove_event (event& ev, event_data* evd_);

void conditionally_remove_event(event& ev);

void update_adjacent_events (event& ev);

void add_reverse_event (event& old_ev);

void perform_event (event& ev);

bool pop_next_event ();

#endif