/** Handling of entities.
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Mar 30, 2020
 *
 * @file
 */

#include <math.h>

#include "global_variables.h"
#include "entity.h"
#include "probability.h"

using namespace std;

/** Add a new entity of a specific type.
 *
 * (the entity is not returned)
 */
void add_entity (
        entity_type et ///< [in] Entity type of the entity to be generated anew.
        )
{
    auto e = ++max_e;
    es.insert(e);
    _e2et[E(e)] = et;
    et2es[et].push_back(e);
    e2label[e] = to_string(e);
    e2outs[e] = { { .rat_out = RT_ID, .e_other = e } };
    e2ins[e] = { { .e_other = e, .rat_in = RT_ID } };
}

/** Return a random entity of a specific type.
 *
 * The entity is drawn uniformly at random.
 */
entity random_entity (
        entity_type et ///< [in] Entity type of which a random entity is needed
        )
{
    int pos = floor(uniform(random_variable) * et2es[et].size());
    auto e = et2es[et][pos];
    return e;
}



