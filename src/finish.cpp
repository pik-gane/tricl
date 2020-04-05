/*
 * finish.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include "debugging.h"
#include "gexf.h"
#include "finish.h"

void finish ()
{
    finish_gexf();
    verify_data_consistency();
}
