// make sure this file is only included once:
#ifndef INC_ENTITY_H
#define INC_ENTITY_H

/** Performance-critical inline functions for handling of entities.
 *
 *  \file
 *
 *  See \ref data_model.h for how a \ref tricl::entity relates to other tricl datatypes.
 */

#include "data_model.h"
#include "probability.h"

using namespace std;

entity add_entity (entity_type et, string label);

/** Return an entity uniformly drawn at random from a specific type.
 *
 *  Mainly used in summary events.
 *
 *  \returns the entity
 */
inline entity random_entity (
        entity_type et  ///< [in] Entity type of which a random entity is needed
        )
{
    int pos = floor(uniform(random_variable) * et2es[et].size());
    auto e = et2es[et][pos];
    return e;
}

#endif



