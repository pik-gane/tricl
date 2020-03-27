/*
 * io.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_IO_H
#define INC_IO_H

#include <ostream>

#include "data_model.h"

using namespace std;

ostream& operator<< (ostream& os, const link_type& lt);

ostream& operator<< (ostream& os, const event_type& evt);

ostream& operator<< (ostream& os, const event& ev);

void dump_data ();

#endif

