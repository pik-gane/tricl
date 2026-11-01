/** Forced events and events at deterministic times (see forced.h).
 *
 *  \file
 */

#include "3rdparty/rapidcsv.h"

#include "global_variables.h"
#include "link.h"
#include "event.h"
#include "actions.h"
#include "forced.h"

string forced_in_filename = "";
bool performing_forced = false;

static vector<pair<timepoint, event>> forced_events;  // sorted by time (stable)
static size_t next_forced = 0;

/** Read the forced events (after initialization, since entities and types must be known). */
void read_forced_events ()
{
    if (forced_in_filename == "") return;
    {
        std::ifstream test(forced_in_filename);
        if (!test.good()) throw "cannot open forced events file \"" + forced_in_filename + "\"";
    }
    rapidcsv::Document doc(forced_in_filename, rapidcsv::LabelParams(0, -1), rapidcsv::SeparatorParams(','));
    vector<string> t_col, ec_col, e1_col, rat_col, e3_col;
    try {
        t_col = doc.GetColumn<string>("t");
        ec_col = doc.GetColumn<string>("event");
        e1_col = doc.GetColumn<string>("source");
        rat_col = doc.GetColumn<string>("relationship");
        e3_col = doc.GetColumn<string>("target");
    } catch (const std::exception& e) {
        throw "forced events file " + forced_in_filename + " must have the columns t, event, source, relationship, target (" + e.what() + ")";
    }
    for (size_t i = 0; i < t_col.size(); i++) {
        string where = " in line " + to_string(i + 2) + " of " + forced_in_filename;
        timepoint t = std::stod(t_col[i]);
        event_class ec;
        if (ec_col[i] == "establish") ec = EC_EST;
        else if (ec_col[i] == "terminate") ec = EC_TERM;
        else throw "forced events must be establish or terminate events" + where;
        if (label2e.count(e1_col[i]) == 0) throw "unknown entity \"" + e1_col[i] + "\"" + where;
        if (label2e.count(e3_col[i]) == 0) throw "unknown entity \"" + e3_col[i] + "\"" + where;
        if (label2rat.count(rat_col[i]) == 0) throw "unknown relationship type \"" + rat_col[i] + "\"" + where;
        entity e1 = label2e[e1_col[i]], e3 = label2e[e3_col[i]];
        auto rat = label2rat[rat_col[i]];
        if ((rat == RT_ID) || r_is_action_type[rat] || (e1 == e3)) throw "forced events must change links of a relationship type between different entities" + where;
        if ((evt_id_of(EC_EST, e2et[e1], rat, e2et[e3]) >= 0) || (evt_id_of(EC_TERM, e2et[e1], rat, e2et[e3]) >= 0))
            throw "forced events must be of link types without dynamics" + where;
        forced_events.push_back({ t, { .ec = ec, .e1 = e1, .rat13 = rat, .e3 = e3 } });
    }
    std::stable_sort(forced_events.begin(), forced_events.end(), [](const pair<timepoint, event>& a, const pair<timepoint, event>& b) { return a.first < b.first; });
    if (!forced_events.empty() && (forced_events[0].first < current_t)) throw "forced events must not happen before the start time";
}

bool any_deterministic_events ()
{
    return any_actions || !forced_events.empty();
}

timepoint next_deterministic_time ()
{
    timepoint tf = (next_forced < forced_events.size()) ? forced_events[next_forced].first : INFINITY;
    timepoint tx = any_actions ? next_expiry() : INFINITY;
    return min(tf, tx);
}

/** Apply a forced event (no log-likelihood, not written to the events output). */
static void apply_forced (event ev)
{
    tricllink l = { ev.e1, ev.rat13, ev.e3 };
    if ((ev.ec == EC_EST) == link_exists(l)) {
        if (!quiet) cout << "skipping forced event " << ev << " at t=" << current_t << " (link " << (link_exists(l) ? "exists already" : "does not exist") << ")" << endl;
        return;
    }
    static event_data forced_evd;
    forced_evd = {};
    forced_evd.effective_rate = 1.0;  // (so that its log-likelihood contribution log(1) is zero)
    performing_forced = true;
    perform_event(ev, &forced_evd);
    performing_forced = false;
}

/** Process all events at deterministic times up to time t (inclusive), in time order
 *  (forced events before expiries at equal times).
 */
void process_deterministic_until (timepoint t)
{
    while (true) {
        timepoint tf = (next_forced < forced_events.size()) ? forced_events[next_forced].first : INFINITY;
        timepoint tx = any_actions ? next_expiry() : INFINITY;
        if (!(min(tf, tx) <= t)) break;
        if (tf <= tx) {
            advance_time(max(current_t, tf));
            apply_forced(forced_events[next_forced++].second);
        } else {
            process_one_expiry();
        }
        if (n_infinite_effective_rates > 0) break;  // immediate events triggered by a forced event happen first
    }
}
