/** Handling of entities.
 *
 *  \file
 *
 *  See \ref data_model.h for how a \ref tricl::entity relates to other tricl datatypes.
 */

#include "global_variables.h"
#include "entity.h"

/** Add a new entity of a specific type.
 *
 *  \returns the new entity
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
    e2ins[e]  = { { .e_source = e, .rat_in = RT_ID } };

    return e;
}

