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

void schedule_event (event& ev, event_data* evd_);

void reschedule_event (event& ev, event_data* evd_);

void add_event (event& ev);

void remove_event (event& ev, event_data* evd_);

void update_adjacent_events (event& ev);

void add_reverse_event (event& old_ev);

void perform_event (event& ev);

bool pop_next_event ();

#endif
