#ifndef INC_ENTITY_H
#define INC_ENTITY_H

#include "data_model.h"

using namespace std;

entity add_entity (
        entity_type et, ///< [in] entity type of the entity to be generated anew.
        string label    ///< [in] entity label (if "", will be generated)
        );

entity random_entity (entity_type et);

#endif



