/*
 * finish.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include "global_variables.h"
#include "debugging.h"
#include "io.h"
#include "gexf.h"
#include "finish.h"

void finish ()
{
    current_t = max_t;
    log_status();
    finish_gexf();
    if (debug) verify_data_consistency();
}
