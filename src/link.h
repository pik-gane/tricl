// make sure this file is only included once:
#ifndef INC_LINK_H
#define INC_LINK_H

#include "data_model.h"

bool link_exists (tricllink& l);

void add_link (tricllink& l);

void delete_link (tricllink& l);

void do_random_link (probability p, entity e1, relationship_or_action_type rat13, entity e3);

#endif
