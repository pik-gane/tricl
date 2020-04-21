// make sure this file is only included once:
#ifndef INC_GEXF_H
#define INC_GEXF_H

#include "data_model.h"

extern unordered_map<tricllink, timepoint> gexf_edge2start;

void init_gexf ();

void gexf_output_edge (tricllink& l);

void finish_gexf ();

#endif
