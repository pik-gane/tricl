/*
 * gexf.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_GEXF_H
#define INC_GEXF_H

#include <fstream>
#include <unordered_map>

#include "data_model.h"

using namespace std;

extern ofstream gexf;
extern unordered_map<link, timepoint> gexf_edge2start;

void init_gexf ();

void gexf_output_edge (link& l);

void finish_gexf ();

#endif
