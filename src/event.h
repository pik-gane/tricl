/*
 * event.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_EVENT_H
#define INC_EVENT_H

#include <assert.h>

#include "data_model.h"

using namespace std;

inline bool event_is_summary(const event& ev) {
    // summary events are encoded by using the negative entity type id as "entity id":
    bool res = ((ev.e1 < 0) || (ev.e3 < 0));
    if (res) {
        assert (ev.e1 < 0);
        assert (ev.e3 < 0);
        assert (ev.ec == EC_EST);
    }
    return res;
}

inline entity_type summary_et1(const event& ev) {
    return (entity_type) -ev.e1;
}
inline entity_type summary_et3(const event& ev) {
    return (entity_type) -ev.e3;
}

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
