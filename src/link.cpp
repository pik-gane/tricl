/** Handling of links
 *
 *  \file
 *
 *  See \ref data_model.h for how a \ref tricl::tricllink relates to other tricl datatypes.
 */

#include "global_variables.h"
#include "probability.h"
#include "event.h"
#include "gexf.h"
#include "link.h"

/** \returns whether link currently exists.
 */
bool link_exists (tricllink& l)
{
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;
    return (e2outs[e1].count({ .rat_out = rat13, .e_target = e3 }) > 0);
}

/** add a link.
 */
void add_link (tricllink& l)
{
    assert ((!link_exists(l)) && (l.e1 != l.e3));
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;

    // keep inleg and outleg sets consistent:
    e2outs[e1].insert({ .rat_out = rat13, .e_target = e3 });
    e2ins[e3].insert({ .e_source = e1, .rat_in = rat13 });

    // register birth time for later output:
    if (rat13 != RT_ID) gexf_edge2start[l] = current_t;

    // update counts:
    lt2n[{e2et[e1], rat13, e2et[e3]}]++;
    n_links++;
}

/** delete a link.
 */
void delete_link (tricllink& l)
{
    assert (link_exists(l));
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;

    // keep inleg and outleg sets consistent:
    e2outs[e1].erase({ .rat_out = rat13, .e_target = e3 });
    e2ins[e3].erase({ .e_source = e1, .rat_in = rat13 });

    // output to gexf:
    if (rat13 != RT_ID) gexf_output_edge(l);

    // update counts:
    lt2n[{e2et[e1], rat13, e2et[e3]}]--;
    n_links--;
}

/** perform an event that generates a random link during installation.
 */
void do_random_link (probability p, entity e1, relationship_or_action_type rat13, entity e3) {
    if (uniform(random_variable) < p) {
        tricllink l = { e1, rat13, e3 };
        if (!link_exists(l)) {
            event ev = { .ec=EC_EST, e1, rat13, e3 };
            conditionally_remove_event(ev);
            // prepair total effective rate since call to perform_event will subtracted INFINITY from it:
            add_effective_rate(INFINITY);
            perform_event(ev, &sure_evd);
        }
    }
}
