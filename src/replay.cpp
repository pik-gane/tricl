/** Replay mode: compute the log-likelihood (and optionally its gradient) of a given sequence of events.
 *
 *  \file
 *
 *  Instead of drawing events at random, the events are read from a csv file with the columns
 *  t, event ("establish" or "terminate"), source, relationship, target (as written by --events-out),
 *  and applied in order. The bookkeeping of event rates is exactly the same as in a simulation,
 *  so the log-likelihood of the given sequence under the model is accumulated as in a simulation:
 *  each event contributes log(its current rate), and each advance of model time contributes
 *  -(total rate) * dt. If limits:t is finite, the log-probability of no further event until then is included.
 *
 *  The event sequence must be complete, i.e. it must also contain the "immediate" events (with infinite rates)
 *  that the model performs in response to other events, and events that are impossible under the model
 *  (e.g. the termination of a non-existing link) are reported as errors.
 *  Companion events of inverse relationship types are implied and must not be listed.
 */

#include "3rdparty/rapidcsv.h"

#include "global_variables.h"
#include "entity.h"
#include "link.h"
#include "event.h"
#include "io.h"
#include "debugging.h"
#include "replay.h"

#define RECOMPUTE_TOTAL_ER_EVERY 65536  ///< No. of events after which the total effective rate is recomputed exactly to remove floating point drift

/** Compose an error message for an event that cannot happen under the model. */
static string impossible (const string& why, size_t line)
{
    return "the event in line " + to_string(line) + " of " + events_in_filename + " is impossible under the model: " + why;
}

/** Replay the events given in the file \ref events_in_filename.
 */
void replay ()
{
    if (!quiet) cout << "REPLAYING events from " << events_in_filename << " ..." << endl;
    {
        std::ifstream test(events_in_filename);
        if (!test.good()) throw "cannot open events file \"" + events_in_filename + "\"";
    }
    rapidcsv::Document doc(events_in_filename,
            rapidcsv::LabelParams(0, -1),  // first row contains column names
            rapidcsv::SeparatorParams(','));
    vector<string> t_col, ec_col, e1_col, rat_col, e3_col;
    try {
        t_col = doc.GetColumn<string>("t");
        ec_col = doc.GetColumn<string>("event");
        e1_col = doc.GetColumn<string>("source");
        rat_col = doc.GetColumn<string>("relationship");
        e3_col = doc.GetColumn<string>("target");
    } catch (const std::exception& e) {
        throw "events file " + events_in_filename + " must have the columns t, event, source, relationship, target (" + e.what() + ")";
    }

    static event_data evd;  // data of the event being performed (must outlive perform_event)
    for (size_t i = 0; i < t_col.size(); i++)
    {
        size_t line = i + 2;  // line number in the file (1-based, after the header)

        // parse the event:
        timepoint t;
        try {
            t = std::stod(t_col[i]);
        } catch (const std::exception&) {
            throw "cannot parse time \"" + t_col[i] + "\" in line " + to_string(line) + " of " + events_in_filename;
        }
        event_class ec;
        if (ec_col[i] == "establish") ec = EC_EST;
        else if (ec_col[i] == "terminate") ec = EC_TERM;
        else throw "unknown event class \"" + ec_col[i] + "\" in line " + to_string(line) + " of " + events_in_filename + " (must be establish or terminate)";
        auto e1_it = label2e.find(e1_col[i]);
        if (e1_it == label2e.end()) throw "unknown entity \"" + e1_col[i] + "\" in line " + to_string(line) + " of " + events_in_filename;
        auto e3_it = label2e.find(e3_col[i]);
        if (e3_it == label2e.end()) throw "unknown entity \"" + e3_col[i] + "\" in line " + to_string(line) + " of " + events_in_filename;
        auto rat_it = label2rat.find(rat_col[i]);
        if (rat_it == label2rat.end()) throw "unknown relationship type \"" + rat_col[i] + "\" in line " + to_string(line) + " of " + events_in_filename;
        entity e1 = e1_it->second, e3 = e3_it->second;
        relationship_or_action_type rat13 = rat_it->second;
        if (rat13 == RT_ID) throw impossible("the identity relationship cannot change", line);
        if (e1 == e3) throw impossible("self-links are not allowed", line);
        event ev = { .ec = ec, .e1 = e1, .rat13 = rat13, .e3 = e3 };

        // advance time:
        if (t < current_t) throw "events in " + events_in_filename + " are not in chronological order (line " + to_string(line) + ")";
        if (n_infinite_effective_rates > 0)
        {
            // an immediate event is pending, so the next event must be one of them and happen right now:
            if (t > current_t) throw impossible("an immediate event is pending at t=" + to_string(current_t) + " but the next event in the file happens later", line);
        }
        else
        {
            advance_time(t);  // accumulates the log-probability that nothing happened in between
        }

        // determine the event's current rate and unregister it, as pop_next_event() would do:
        int evt_id = evt_id_of(ec, e2et[e1], rat13, e2et[e3]);
        auto it = ev2data.find(ev);
        if (it != ev2data.end())  // event is scheduled individually
        {
            evd = it->second;
            if (!(evd.effective_rate > 0)) throw impossible("its current rate is zero", line);
            if ((n_infinite_effective_rates > 0) && (evd.effective_rate < INFINITY)) throw impossible("an immediate event is pending", line);
            add_event_log_gradient(&evd, evt_id);
            remove_event(ev, &it->second);
        }
        else if (ec == EC_EST)  // event may be covered by a summary event
        {
            tricllink l = { e1, rat13, e3 };
            if (link_exists(l)) throw impossible("the link exists already", line);
            if ((evt_id < 0) || !(evtid2summary_single_er[evt_id] > 0)) throw impossible("this type of link cannot be established spontaneously", line);
            if (n_infinite_effective_rates > 0) throw impossible("an immediate event is pending", line);
            evd = {};
            evd.effective_rate = evtid2summary_single_er[evt_id];  // (legs cannot influence establishment yet)
            evd.t = current_t;
            add_summary_event_log_gradient(evt_id);
            subtract_summary_shares(evt_id, 1);
        }
        else
        {
            throw impossible("the link does not exist", line);
        }

        // perform it:
        current_ev = ev;
        current_evd_ = &evd;
        ++n_events;
        log_state();
        perform_event(ev, &evd);
        if (n_events % RECOMPUTE_TOTAL_ER_EVERY == 0) total_finite_effective_rate = max(0.0, compute_total_finite_er());
        if (n_events >= max_n_events) break;
    }

    // the log-probability that nothing happens until the time limit
    // (as in a simulation, the observation ends with the last event if the event limit was reached):
    if ((max_t < INFINITY) && (n_events < max_n_events))
    {
        if (max_t < current_t) throw "the time limit limits:t lies before the last event in " + events_in_filename;
        finish_time();
    }
    if (!quiet) cout << "...REPLAY FINISHED after " << n_events << " events." << endl;
}
