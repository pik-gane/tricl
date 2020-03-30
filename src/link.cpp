/*
 * link.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include <assert.h>

#include "global_variables.h"
#include "probability.h"

#include "event.h"

#include "gexf.h"

#include "link.h"

using namespace std;

bool link_exists (link& l)
{
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;
    return (e2outs[e1].count({ .e = e3, .r = rat13}) > 0);
}

void add_link (link& l)
{
    assert ((!link_exists(l)) && ((l.e1 != l.e3) || (l.rat13 == RT_ID)));
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;
    e2outs[e1].insert({ .e = e3, .r = rat13 });
    e2ins[e3].insert({ .e = e1, .r = rat13 });
    if (rat13 != RT_ID) gexf_edge2start[l] = current_t;
    n_links++;
}

void del_link (link& l)
{
    assert (link_exists(l));
    auto e1 = l.e1, e3 = l.e3; auto rat13 = l.rat13;
    e2outs[e1].erase({ .e = e3, .r = rat13 });
    e2ins[e3].erase({ .e = e1, .r = rat13 });
    if (rat13 != RT_ID) gexf_output_edge(l);
    n_links--;
}

void do_random_link (probability p, entity_type e1, relationship_or_action_type rat13, entity_type e3) {
    if (uniform(random_variable) < p) {
        link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 };
        if (!link_exists(l)) {
            event ev = { .ec = EC_EST, .e1 = e1, .rat13 = rat13, .e3 = e3 };
            auto evd_ = &ev2data[ev];
            if (evd_->t > -INFINITY) remove_event(ev, evd_);
            perform_event(ev);
        }
    }
}

