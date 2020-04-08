/*
 * io.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include <iostream>
#include <ostream>
#include "3rdparty/rapidcsv.h"

#include "data_model.h"
#include "global_variables.h"
#include "entity.h"
#include "io.h"

using namespace std;

// I/O:

ostream& operator<< (ostream& os, const link_type& lt) {
    os << et2label[lt.et1]
        << " " << rat2label[lt.rat13] << " "
        << et2label[lt.et3];
    return os;
}

ostream& operator<< (ostream& os, const event_type& evt) {
    os << ec2label[evt.ec]
            << " \"" << et2label[evt.et1]
            << " " << rat2label[evt.rat13] << " "
            << et2label[evt.et3]
            << "\"";
    return os;
}

ostream& operator<< (ostream& os, const event& ev) {
    if (ev.e1 < 0) {
        os << ec2label[ev.ec]
                << " \"some " << et2label[-ev.e1]
                << " " << rat2label[ev.rat13] << " some "
                << et2label[-ev.e3]
                << "\"";
    } else {
        os << ec2label[ev.ec]
                << " \"" << e2label[ev.e1]
                << " " << rat2label[ev.rat13] << " "
                << e2label[ev.e3]
                << "\"";
    }
    return os;
}

ostream& operator<< (ostream& os, const event_data& evd) {
    os << "na=" << evd.n_angles << " ar=" << evd.attempt_rate << " pu=" << evd.success_probunits << " t=" << evd.t;
    return os;
}

void read_links_csv (
        string filename,    ///< name of infile
        int skip_rows,      ///< no. of rows to skip at start of file (e.g. header lines)
        int max_rows,       ///< no. of rows to read at most
        char delimiter,     ///< delimiter, e.g. " " or "," or "\t"
        int e1_col,         ///< no. of column containing e1 labels
        int rat13_col,      ///< no. of column containing rat13 labels (or -1 if type is fixed to value of "rat13_fixed")
        int e3_col,         ///< no. of column containing e3 labels
        entity_type et1_default, ///< default entity-type for previously unregistered e1s (if -1, no new e1s are allowed)
        relationship_or_action_type rat13_fixed, ///< fixed type if rat13_col == -1 (otherwise NO_RAT)
        entity_type et3_default, ///< default entity-type for previously unregistered e3s (if -1, no new e3s are allowed)
        string e_prefix     ///< prefix to prepend to entity labels before storing them
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
        auto e1label = e_prefix + e1labels[row], e3label = e_prefix + e3labels[row];
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

void dump_data () // dump important data to stdout for debugging
{
    dump_links();
    cout << "t2be:" << endl;
    for (auto& [t, ev] : t2be) cout << " " << t << ": " << ev << endl;
    cout << "ev2data:" << endl;
    for (auto& [ev2, evd2] : ev2data) cout << " " << ev2 << ": " << evd2 << endl;
}


