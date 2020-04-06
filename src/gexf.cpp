/*
 * gexf.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include "data_model.h"
#include "global_variables.h"
#include "gexf.h"

using namespace std;

ofstream gexf;
unordered_map<link, timepoint> gexf_edge2start = {};

void init_gexf () {
    gexf.open(gexf_filename);
    gexf << R"V0G0N(<?xml version="1.0" encoding="UTF-8"?>
<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">
    <meta>
        <creator>tricl</creator>
        <description>dynamic graph generated by tricl model</description>
    </meta>
    <graph mode="dynamic" defaultedgetype="directed">
        <attributes class="node">
            <attribute id="0" title="entity type" type="string"/>
        </attributes>
        <attributes class="edge">
            <attribute id="1" title="relationship_or_action_type" type="string"/>
        </attributes>
        <nodes>
)V0G0N";
    for (auto& e : es) {
        gexf << "<node id=\"" << e << "\" label=\"" << e2label[e]
             << "\" start=\"0.0\" end=\"" << max_t
             << "\"><attvalues><attvalue for=\"0\" value=\""
             << et2label[_e2et[E(e)]] << "\"/></attvalues></node>" ;
    }
    gexf << R"V0G0N(
        </nodes>
        <edges>
)V0G0N";
}

void gexf_output_edge (link& l) {
    entity e1 = l.e1, e3 = l.e3; relationship_or_action_type rat13 = l.rat13;
    if (rat13 != RT_ID) {
        gexf << "\t\t\t<edge id=\"" << e1 << "_" << rat13 << "_" << e3 << "_" << current_t
             << "\" source=\"" << e1 << "\" target=\"" << e3
             << "\" start=\"" << gexf_edge2start.at(l) << "\" end=\"" << current_t
             << "\"><attvalues><attvalue for=\"1\" value=\"" << rat2label[rat13]
             << "\"/></attvalues></edge>" << endl;
    }
    gexf_edge2start.erase(l);
}

void finish_gexf () {
    current_t = max_t;
    for (auto& [e1, outs] : e2outs) {
        for (auto& [rat13, e3] : outs) {
            link l = { .e1 = e1, .rat13 = rat13, .e3 = e3 };
            gexf_output_edge(l);
        }
    }
    gexf << R"V0G0N(
        </edges>
    </graph>
</gexf>
)V0G0N";
    gexf.close();
}
