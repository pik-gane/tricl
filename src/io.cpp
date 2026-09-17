/** Handling of input/output other than gexf and graphviz.
 *
 *  \file
 */

// use header-only csv library:
#include "3rdparty/rapidcsv.h"

#include <iomanip>
#include <chrono>

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
    os << "na=" << evd.n_angles << " ar=" << total_attempt_rate(&evd) << " spu=" << total_success_probunits(&evd) << " er=" << evd.effective_rate << " t=" << evd.t;
    return os;
}

} // end of namespace tricl


int n_er_times_dt = 0;
double sum_er_times_dt = 0, sum_er_times_dt_squared = 0;
double avg_er_times_dt = 0;

/** Output simple statistics for current model state to stdout.
 */
void log_state (
        bool final  ///< [in] whether this is the final state of the simulation
        )
{
    if (silent) return;
    if (debug) if (last_dt > 0) {  // TODO find out why this product is not 1 on average!
        n_er_times_dt++;
        double f = total_finite_effective_rate * last_dt;
        sum_er_times_dt += f;
        sum_er_times_dt_squared += f*f;
        double avg = sum_er_times_dt / n_er_times_dt,
                stddev = sqrt(sum_er_times_dt_squared / n_er_times_dt - avg*avg),
                stderr = stddev / sqrt(n_er_times_dt);
        avg_er_times_dt = 0.99 * avg_er_times_dt + 0.01 * f;
        cout << "!! " << avg << " +- " << stderr << ", " << avg_er_times_dt << ", er: " << total_finite_effective_rate << " =? " << compute_total_finite_er() << endl;
    }
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
        // in quiet mode, the status line is only updated a few times per second (formatting it for every event is costly):
        static auto last_output_time = std::chrono::steady_clock::now() - std::chrono::hours(1);
        auto now = std::chrono::steady_clock::now();
        if ((!final) && (now - last_output_time < std::chrono::milliseconds(200))) return;
        last_output_time = now;
        cout << fixed << n_events << ": logl " << cumulative_logl << ", er " << total_finite_effective_rate << ", ld " << ld << ", ad " << ad << ", q " << q << ".  t " << current_t << (final ? "\n" : "\r") << std::flush;
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
        cout << " | stats: logl " << cumulative_logl << ", er " << total_finite_effective_rate << ", ld " << ld << ", ad " << ad << ", q " << q << endl;
        if (current_t < max_t) cout << "at t=" << current_t << " " << current_ev << defaultfloat << endl;
    }
    else
    {
        cout << fixed << n_events << ": logl " << cumulative_logl << ", er " << total_finite_effective_rate << ", ld " << ld << ", ad " << ad << ", q " << q;
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

        int e1_col,         ///< [in] no. of column containing e1 entity labels
        string e1_prefix,   ///< [in] prefix to prepend to e1 labels before storing them
        int et1_col,        ///< [in] no. of column containing et1 entity type labels (or -1 if et1_default is used instead)
        entity_type et1_default, ///< [in] default entity-type for previously unregistered e1s (if -1, no new e1s are allowed or et1_col is used instead)

        int rat13_col,      ///< [in] no. of column containing rat13 labels (or -1 if type is fixed to value of "rat13_fixed")
        relationship_or_action_type rat13_fixed, ///< [in] fixed type if rat13_col == -1 (otherwise NO_RAT)

        int e3_col,         ///< [in] no. of column containing e3 entity labels
        string e3_prefix,    ///< [in] prefix to prepend to e3 labels before storing them
        int et3_col,        ///< [in] no. of column containing et3 entity type labels (or -1 if et3_default is used instead)
        entity_type et3_default  ///< [in] default entity-type for previously unregistered e3s (if NO_ET, no new e3s are allowed or et3_col is used instead)
        )
{
    {
        std::ifstream test(filename);
        if (!test.good()) throw "cannot open initial links file \"" + filename + "\"";
    }
    rapidcsv::Document doc(filename,
            rapidcsv::LabelParams(-1, -1), // access all rows and cols
            rapidcsv::SeparatorParams(delimiter)
        );

    auto e1labels = doc.GetColumn<string>(e1_col);
    auto e3labels = doc.GetColumn<string>(e3_col);
    int nrows = e1labels.size();
    vector<string> et1labels, rat13labels, et3labels;
    if (et1_col > -1) et1labels = doc.GetColumn<string>(et1_col);
    if (rat13_col > -1) rat13labels = doc.GetColumn<string>(rat13_col);
    if (et3_col > -1) et3labels = doc.GetColumn<string>(et3_col);
    for (int row = skip_rows; row < min(nrows, skip_rows+max_rows); row++) {
        if (debug) cout << "row " << row << " " << e1labels[row] << " " << e3labels[row] << endl;
        auto e1label = e1_prefix + e1labels[row], e3label = e3_prefix + e3labels[row];
        entity_type et1, et3;
        entity e1, e3;

        // get and verify et1:
        if (label2e.count(e1label) == 0)
        {
            if ((et1_col > -1) && (et1labels[row] != "")) {
                if (label2et.count(et1labels[row])) et1 = label2et.at(et1labels[row]);
                else throw "unknown et1 label " + et1labels[row];
            }
            else if (et1_default != NO_ET) et1 = et1_default;
            else throw "cannot determine entity type for " + e1label;
            e1 = add_entity(et1, e1label);
            if (verbose) cout << "  entity " << e1 << ": " << et2label.at(et1) << ": " << e1label << endl;
        }
        else
        {
            e1 = label2e.at(e1label);
            et1 = e2et[e1];
            if ((et1_col > -1) && (et1labels[row] != "") && (et1labels[row] != et2label[et1])) throw
                    "wrong et1 label " + et1labels[row] + " for known entity " + e1label;
        }

        // get and verify et3:
        if (label2e.count(e3label) == 0)
        {
            if ((et3_col > -1) && (et3labels[row] != "")) {
                if (label2et.count(et3labels[row])) et3 = label2et.at(et3labels[row]);
                else throw "unknown et3 label " + et3labels[row];
            }
            else if (et3_default != NO_ET) et3 = et3_default;
            else throw "cannot determine entity type for " + e3label;
            e3 = add_entity(et3, e3label);
            if (verbose) cout << "  entity " << e3 << ": " << et2label.at(et3) << ": " << e3label << endl;
        }
        else
        {
            e3 = label2e.at(e3label);
            et3 = e2et[e3];
            if ((et3_col > -1) && (et3labels[row] != "") && (et3labels[row] != et2label[et3])) throw
                    "wrong et3 label " + et3labels[row] + " for known entity " + e3label;
        }
        if (debug) cout << e1 << " " << e3 << endl;
        auto rat13 = (rat13_col >= 0) ? label2rat.at(rat13labels[row]) : rat13_fixed;
        if (debug) cout << rat13 << endl;
        initial_links.insert({ e1, rat13, e3 });
    }
}

// output of performed events to a csv file:

ofstream events_out;  ///< stream for the optional csv output of all performed events

/** Quote a string for csv output if necessary.
 */
static string csv_quote (const string& s)
{
    if (s.find_first_of(",\"\n") == string::npos) return s;
    string res = "\"";
    for (char c : s) {
        if (c == '"') res += "\"\""; else res += c;
    }
    return res + "\"";
}

/** Open the csv file for performed events (if requested) and write its header.
 *
 *  Columns: model time, event class ("establish" or "terminate"), source entity label,
 *  relationship type label, target entity label.
 *  Companion events for inverse relationship types are not written since they are implied.
 */
void open_events_out ()
{
    if (events_out_filename == "") return;
    events_out.open(events_out_filename);
    if (!events_out.good()) throw "cannot open events output file \"" + events_out_filename + "\"";
    events_out << std::setprecision(17) << "t,event,source,relationship,target" << endl;
}

/** Write a performed event to the csv file (if open).
 */
void write_event_out (const event& ev)
{
    if (!events_out.is_open()) return;
    events_out << current_t
            << "," << ((ev.ec == EC_EST) ? "establish" : (ev.ec == EC_TERM) ? "terminate" : "act")
            << "," << csv_quote(e2label[ev.e1])
            << "," << csv_quote(rat2label[ev.rat13])
            << "," << csv_quote(e2label[ev.e3]) << "\n";
}

/** Close the csv file for performed events (if open).
 */
void close_events_out ()
{
    if (events_out.is_open()) events_out.close();
}

/** Output a one-line JSON summary of the final state to stdout.
 */
void output_json_summary ()
{
    cout << std::setprecision(17)
         << "{\"events\": " << n_events
         << ", \"t\": " << current_t
         << ", \"logl\": " << cumulative_logl
         << ", \"links\": " << n_links
         << ", \"angles\": " << n_angles
         << ", \"total_rate\": " << total_finite_effective_rate
         << ", \"seed\": " << seed
         << "}" << endl;
}

/** (for debugging purposes)
 */
void dump_links () {
    cout << "e2outs:" << endl;
    for (entity e1 = 1; e1 <= max_e; e1++) {
        for (auto& l : e2outs[e1]) cout << " " << e2label[e1] << " " << rat2label[l.rat_out] << " " << e2label[l.e_target] << endl;
    }
    cout << "e2ins:" << endl;
    for (entity e3 = 1; e3 <= max_e; e3++) {
        for (auto& l : e2ins[e3]) cout << " " << e2label[l.e_source] << " " << rat2label[l.rat_in] << " " << e2label[e3] << endl;
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

