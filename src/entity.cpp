/*
 * entity.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */


#include <math.h>

#include "global_variables.h"
#include "entity.h"
#include "probability.h"

using namespace std;

void add_entity (entity_type et)
{
    auto e = ++max_e;
    es.insert(e);
    _e2et[E(e)] = et;
    et2es[et].push_back(e);
    e2label[e] = to_string(e);
    e2outs[e] = e2ins[e] = { { .e = e, .r = RT_ID } };
}
entity random_entity (entity_type et)
{
    int pos = floor(uniform(random_variable) * et2es[et].size());
    auto e = et2es[et][pos];
    return e;
}



