/*
 * link.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_LINK_H
#define INC_LINK_H

#include "data_model.h"

bool link_exists (link& l);

void add_link (link& l);

void del_link (link& l);

void do_random_link (probability p, entity_type e1, relationship_or_action_type rat13, entity_type e3);

#endif
