/** Generate some overview diagrams using graphviz
 *
 *  \file
 *
 *  This depends on <https://www.graphviz.org/> and will generate pdf files.
 *  Filenames start with the prefix specified in the config file
 *  as "diagram prefix".
 */

/*
 * TODO:
 * - action types dashed
 */

#include "global_variables.h"
#include "graphviz.h"

ofstream dot;  ///< the respective output streams

/** Draw a type diagram containing all entity and relationship or action types
 */
void do_graphviz_type_diagram ()
{
    auto dotname = diagram_fileprefix + "_types.dot",
            pdfname = diagram_fileprefix + "_types.pdf";

    // compose graphviz input file in "dot" language:
    dot.open(dotname);
    dot << "digraph types {" << endl;
    dot << " rankdir=\"LR\"" << endl << " ranksep=2.0" << endl;
    for (auto& [et, label] : et2label) {
        dot << " " << et << " [label=\"" << label << "\"]" << endl;
    }
    for (auto& [lt, dummy] : lt2n) {
        if (rat2inv.at(lt.rat13) != lt.rat13) { // non-symmetric
            dot << " " << lt.et1 << " -> " << lt.et3 << " [label=\"" << rat2label.at(lt.rat13) << "\"]" << endl;
        } else if (lt.et1 <= lt.et3) { // symmetric, draw only one edge
            dot << " " << lt.et1 << " -> " << lt.et3 << " [label=\"" << rat2label.at(lt.rat13) << "\", arrowhead=none]" << endl;
        }
    }
    dot << "}" << endl;
    dot.close();

    // use graphviz program "dot" to render a pdf:
    string prg = "dot", cmd;
    if (true || debug) {
        cmd = prg + " -v -Tpdf -o\"" + pdfname + "\" \"" + dotname + "\" >\"" + diagram_fileprefix + "_types.log\" 2>&1";
        if (!quiet) cout << " rendering type diagram via command: " << cmd << endl;
    } else {
        cmd = prg + " -Tpdf -o\"" + pdfname + "\" \"" + dotname + "\"";
    }
    if (system(cmd.c_str()) != 0)
        cout << "WARNING: could not render type diagram. Is graphviz installed?" << endl;;
}

/** Draw a modularized diagram showing the dynamic rules of the model.
 *
 *  Each link type gets one panel
 *  with two subpanels for establishment and termination events for this link type.
 *  Each subpanel shows all angles influencing the event,
 *  annotated with the attempt rate and success probability units
 *  contributed by this influence type.
 */
void do_graphviz_dynamics_diagram ()
{
    auto dotname = diagram_fileprefix + "_dynamics.dot",
            pdfname = diagram_fileprefix + "_dynamics.pdf";

    // compose graphviz input file in "dot" language:
    dot.open(dotname);
    dot << "digraph dynamics {" << endl;
    dot << " rankdir=\"LR\"" << endl << " ranksep=0.5" << endl;
    for (auto& [lt, dummy] : lt2n) {
        auto et1 = lt.et1, et3 = lt.et3;
        auto rat13 = lt.rat13;
        string sub = "cluster_" + to_string(et1) + "_" + to_string(rat13) + "_" + to_string(et3),
                et1l = et2label.at(et1), et3l = et2label.at(et3),
                rat13l = rat2label.at(rat13);
        dot << " subgraph " << sub << " { label=\"" << et1l << " " << rat13l << " " << et3l << "\"; fontsize=21.0; fontname=\"times bold\"; style=filled; color=lightgrey;" << endl;
        for (auto& ec : vector<event_class> {EC_EST, EC_TERM}) {
            event_type evt = {ec, et1, rat13, et3};
            if (evt2base_probunits.count(evt) == 0) continue;
            string subsub = sub + "_" + to_string(ec),
                    v1 = subsub + "_source", v3 = subsub + "_target";
            dot << "  subgraph " << subsub << " { label=\"" << ((ec == EC_EST) ? "establish" : "terminate") << "\"; style=filled; color=white; fontsize=17.0; fontname=\"times bold italic\"" << endl;
            dot << "   " << v1 << " [label=\"" << et2label.at(et1) << "\"style=filled; color=lightgrey]" << endl;
            dot << "   " << v3 << " [label=\"" << et2label.at(et3) << "\"style=filled; color=lightgrey]" << endl;
            // in dynamics diagrams, also symmetric links have an arrow head:
            dot << "   " << v1 << " -> " << v3 << " [label=\"" << rat13l << "\\n";
            auto ar = (evt2base_attempt_rate.count(evt) > 0) ? evt2base_attempt_rate.at(evt) : 0.0,
                    spu = evt2base_probunits.at(evt),
                    til = evt2left_tail.at(evt), tir = evt2right_tail.at(evt);
            if (ar > 0.0) dot << "AR +" << ar << ", ";
            dot << "TI " << til;
            if (tir != til) dot << "/" << tir;
            if (spu != 0.0) dot << ", SPU " << ((spu > 0.0) ? "+" : "") << spu;
            dot << "\"; style=bold]" << endl;
            for (auto& [inflt, ar] : inflt2attempt_rate) {
                if (inflt.evt == evt) {
                    auto rat12 = inflt.at.rat12, rat23 = inflt.at.rat23;
                    auto et2 = inflt.at.et2;
                    string v2 = subsub + "_" + to_string(rat12) + "_" + to_string(et2) + "_" + to_string(rat23);
                    dot << "   " << v2 << " [label=\"" << et2label.at(et2) << "\\nAR +" << ar;
                    if (inflt2delta_probunits.count(inflt) > 0) {
                        auto dspu = inflt2delta_probunits.at(inflt);
                        if (dspu != 0.0) dot << ", " << "SPU" << ((dspu > 0.0) ? "+" : "") <<  dspu;
                    }
                    dot << "\"]" << endl;
                    dot << "   " << v1 << " -> " << v2 << " [label=\"" << rat2label.at(rat12) << "\"]";
                    dot << "   " << v2 << " -> " << v3 << " [label=\"" << rat2label.at(rat23) << "\"]";
                }
            }
            for (auto& [inflt, dspu] : inflt2delta_probunits) {
                if ((inflt.evt == evt) && (inflt2attempt_rate.count(inflt) == 0)) {
                    auto rat12 = inflt.at.rat12, rat23 = inflt.at.rat23;
                    auto et2 = inflt.at.et2;
                    string v2 = subsub + "_" + to_string(rat12) + "_" + to_string(et2) + "_" + to_string(rat23);
                    dot << "   " << v2 << " [label=\"" << et2label.at(et2);
                    auto dspu = inflt2delta_probunits.at(inflt);
                    if (dspu != 0.0) dot << "\\nSPU " << ((dspu > 0.0) ? "+" : "") <<  dspu;
                    dot << "\"]" << endl;
                    dot << "   " << v1 << " -> " << v2 << " [label=\"" << rat2label.at(rat12) << "\"]";
                    dot << "   " << v2 << " -> " << v3 << " [label=\"" << rat2label.at(rat23) << "\"]";
                }
            }
            dot << "  }" << endl;
        }
        dot << " }" << endl;
    }
    // footnote legend:
    dot << " label = \"* AR = attempt rate, TI = tail indices, SPU = success probability units\"; labelloc = \"b\";" << endl << "}" << endl;
    dot.close();

    // use graphviz program "dot" to render a pdf:
    string prg = "dot", cmd;
    if (true || debug) {
        cmd = prg + " -v -Tpdf -o\"" + pdfname + "\" \"" + dotname + "\" >\"" + diagram_fileprefix + "_dynamics.log\" 2>&1";
        if (!quiet) cout << " rendering dynamics diagram via command: " << cmd << endl;
    } else {
        cmd = prg + " -Tpdf -o\"" + pdfname + "\" \"" + dotname + "\"";
    }
    if (system(cmd.c_str()) != 0)
        cout << "WARNING: could not render dynamics diagram. Is graphviz installed?" << endl;;
}

/** Generate all diagrams.
 */
void do_graphviz_diagrams ()
{
    if (diagram_fileprefix != "") {
        do_graphviz_type_diagram();
        do_graphviz_dynamics_diagram();
    }
}
