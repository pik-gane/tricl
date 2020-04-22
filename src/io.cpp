/** Handling of input/output other than gexf and graphviz.
 *
 *  \file
 */

// use header-only csv library:
#include "3rdparty/rapidcsv.h"

#include "global_variables.h"
#include "entity.h"
#include "event.h"
#include "io.h"

// overloaded streaming operators need to live in our special namespace:
namespace tricl {

ostream& operator<< (ostream& os, const link_type& lt) {
    os << et2label[lt.et1]
        << " " << rat2label[lt.rat13] << " "
        << et2label[lt.et3];
    return os;
}

ostream& operator<< (ostream& os, const event_type& evt)
{
    os << ec2label[evt.ec]
            << " " << et2label[evt.et1]
            << " " << rat2label[evt.rat13] << " "
            << et2label[evt.et3];
    return os;
}

ostream& operator<< (ostream& os, const angle_type& at) {
    os << rat2label[at.rat12] << " " << et2label[at.et2] << " " << rat2label[at.rat23];
    return os;
}

ostream& operator<< (ostream& os, const event& ev) {
    if (event_is_summary(ev)) {
        os << ec2label[ev.ec]
                << " some " << et2label[summary_et1(ev)]
                << " " << rat2label[ev.rat13] << " some "
                << et2label[summary_et3(ev)];
    } else {
        os << ec2label[ev.ec]
                << " " << e2label[ev.e1]
                << " " << rat2label[ev.rat13] << " "
                << e2label[ev.e3];
    }
    return os;
}

ostream& operator<< (ostream& os, const event_data& evd) {
    os << "na=" << evd.n_angles << " ar=" << evd.attempt_rate << " pu=" << evd.success_probunits << " t=" << evd.t;
    return os;
}

} // end of namespace tricl


/** Output simple statistics for current model state to stdout.
 */
void log_state ()
{
    double
        ne = (double) max_e,                 ///< no. past events
        nl = (double) n_links,               ///< current no. (non-id.) links
        na = (double) n_angles,              ///< current no. angles influencing at least one event
        ld = nl / (ne * ne),                 ///< overall link density
        ad = na / (ne * ne * ne),            ///< overall angle density
        q = (ld > 0.0) ? ad / (ld*ld) : 0.0  ///< quotient between angle density and "expected angle density" (square of link density)
        ;
    if (quiet)
    {
        cout << fixed << n_events << ": ld " << ld << ", ad " << ad << ", q " << q << ".  t " << current_t << "\r";
    }
    else if (lt2n.size() > 1)
    {
        cout << endl << fixed << n_events;
        for (auto& [lt, n] : lt2n) {
            auto rat13 = lt.rat13, rat31 = rat2inv.at(rat13);
            if ((n > 0) && ((rat31 == NO_RAT) || (rat31 >= rat13))) {
                cout << " | " << n << " " << lt;
            }
        }
        cout << " | stats: ld " << ld << ", ad " << ad << ", q " << q << endl;
        if (current_t < max_t) cout << "at t=" << current_t << " " << current_ev << defaultfloat << endl;
    }
    else
    {
        cout << fixed << n_events << ": ld " << ld << ", ad " << ad << ", q " << q;
        if (current_t < max_t) cout << ".  t " << current_t << ": " << current_ev << defaultfloat;
        cout << endl;
    }
}

/** Read initial links from a csv file.
 *
 *  The file could either contain links of just one type,
 *  in which case you specify two columns (containing the source and target entity labels);
 *  or links of various types,
 *  in which case you also specify a column for the relationship or action type.
 *
 *  If entity labels occur that have not been registered, they are assigned to
 *  the entity type given by et1_default or et3_default, respectively.
 *
 *  Entity labels from the file are prefixed by e_prefix before storing them in tricl,
 *  e.g. the file might contain label "1" and you want it to be treated as "agent 1",
 *  so put e1_prefix = "agent ";
 */
void read_links_csv (
        string filename,    ///< [in] name of infile
        int skip_rows,      ///< [in] no. of rows to skip at start of file (e.g. header lines)
        int max_rows,       ///< [in] no. of rows to read at most
        char delimiter,     ///< [in] delimiter, e.g. " " or "," or "\t"
        int e1_col,         ///< [in] no. of column containing e1 labels
        int rat13_col,      ///< [in] no. of column containing rat13 labels (or -1 if type is fixed to value of "rat13_fixed")
        int e3_col,         ///< [in] no. of column containing e3 labels
        entity_type et1_default, ///< [in] default entity-type for previously unregistered e1s (if -1, no new e1s are allowed)
        relationship_or_action_type rat13_fixed, ///< [in] fixed type if rat13_col == -1 (otherwise NO_RAT)
        entity_type et3_default, ///< [in] default entity-type for previously unregistered e3s (if -1, no new e3s are allowed)
        string e1_prefix,   ///< [in] prefix to prepend to e1 labels before storing them
        string e3_prefix    ///< [in] prefix to prepend to e3 labels before storing them
        )
{
    rapidcsv::Document doc(filename,
            rapidcsv::LabelParams(-1, -1), // access all rows and cols
            rapidcsv::SeparatorParams(delimiter)
        );

    auto e1labels = doc.GetColumn<string>(e1_col);
    auto e3labels = doc.GetColumn<string>(e3_col);
    int nrows = e1labels.size();
    std::vector<string> rat13labels;
    if (rat13_col >= 0) rat13labels = doc.GetColumn<string>(rat13_col);
    for (int row = skip_rows; row < min(nrows, skip_rows+max_rows); row++) {
        if (debug) cout << "row " << row << " " << e1labels[row] << " " << e3labels[row] << endl;
        auto e1label = e1_prefix + e1labels[row], e3label = e3_prefix + e3labels[row];
        entity e1, e3;
        if (label2e.count(e1label) == 0) {
            if (et1_default == -1) throw "unknown entity " + e1label;
            e1 = add_entity(et1_default, e1label);
            if (verbose) cout << "  entity " << e1 << ": " << et2label.at(et1_default) << ": " << e1label << endl;
        } else {
            e1 = label2e.at(e1label);
        }
        if (label2e.count(e3label) == 0) {
            if (et3_default == -1) throw "unknown entity " + e3label;
            e3 = add_entity(et3_default, e3label);
            if (verbose) cout << "  entity " << e3 << ": " << et2label.at(et3_default) << ": " << e3label << endl;
        } else {
            e3 = label2e.at(e3label);
        }
        if (debug) cout << e1 << " " << e3 << endl;
        auto rat13 = (rat13_col >= 0) ? label2rat.at(rat13labels[row]) : rat13_fixed;
        if (debug) cout << rat13 << endl;
        initial_links.insert({ e1, rat13, e3 });
    }
}

/** (for debugging purposes)
 */
void dump_links () {
    cout << "e2outs:" << endl;
    for (auto& [e1, outs1] : e2outs) {
        for (auto& [rat13, e3] : outs1) cout << " " << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << endl;
    }
    cout << "e2ins:" << endl;
    for (auto& [e3, ins3] : e2ins) {
        for (auto& [e1, rat13] : ins3) cout << " " << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << endl;
    }
}

/** (for debugging purposes)
 */
void dump_data ()
{
    dump_links();
    cout << "t2be:" << endl;
    for (auto& [t, ev] : t2ev) cout << " " << t << ": " << ev << endl;
    cout << "ev2data:" << endl;
    for (auto& [ev2, evd2] : ev2data) cout << " " << ev2 << ": " << evd2 << endl;
}

