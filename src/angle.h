#ifndef INC_ANGLE_H
#define INC_ANGLE_H

#include "data_model.h"

void add_or_delete_angle (
        event_class ec,
        entity e1, entity_type et1, relationship_or_action_type rat12,
        entity e2, entity_type et2, relationship_or_action_type rat23,
        entity e3, entity_type et3);

angles leg_intersection (entity e1, outleg_set& out1, inleg_set& in3, entity e3);

#endif
