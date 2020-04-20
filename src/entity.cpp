/** Handling of entities.
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Mar 30, 2020
 *
 * @file
 */

#include <iostream>
#include <math.h>

#include "global_variables.h"
#include "entity.h"
#include "probability.h"

/** Add a new entity of a specific type.
 *
 * (the entity is not returned)
 */
entity add_entity (
        entity_type et, ///< [in] entity type of the entity to be generated anew.
        string elabel   ///< [in] entity label (if "", will be generated)
        )
{
    // generate id and register in list:
    auto e = ++max_e;
    es.insert(e);

    // register type:
    e2et[e] = et;
    et2es[et].push_back(e);

    // register label:
    if (elabel == "") {
        elabel = et2label.at(et) + " " + to_string(e);
    }
    if (label2e.count(elabel) > 0) throw "entity label \"" + elabel + "\" is already in use";
    e2label[e] = elabel;
    label2e[elabel] = e;

    // register identity relation:
    e2outs[e] = { { .rat_out = RT_ID, .e_target = e } };
    e2ins[e] = { { .e_source = e, .rat_in = RT_ID } };

    return e;
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



