/** An efficient data structure for storing the schedule.
 *
 *  \file
 */


#include "global_variables.h"


/** An efficient container for pairs of (event, event_data) that supports efficient
 * - lookup, deletion, insertion, and update of event_data by event
 * - lookup of item with minimal event_data.t
 */
class schedule_t
{
    vector<pair<event, timepoint>> ev_t_now = {};  ///< the small linear container for schedule class SC_NOW
    int n_ev_now = 0;

    map<timepoint, event> t2ev_sooner = {};  ///< the tree-type container for schedule class SC_TREE

    unordered_map<event, timepoint> ev2t_later = {};  ///< the hash table for schedule class SC_HASH
    timepoint min_t_later = INFINITY;                 ///< the smallest t of schedule class S_LIST

    unordered_map<event, event_data> ev2data = {};  ///< the hash table for all events' data

    inline auto _find_now (const event& ev)
    {
        auto it = ev_t_now.begin();
        for (; it != ev_t_now.end(); it++)
        {
            if (it->first == ev) return it;
        }
        return it;
    }
    inline event_data at (const event& ev)
    {
        return ev2data.at(ev);
    }
    inline void _update_min_t_later(timepoint removed_t)
    {
        if (removed_t == min_t_later) {
            min_t_later = INFINITY;
            for (const auto& [ev, t] : ev2t_later) if (t < min_t_later) min_t_later = t;
        }
    }
    inline event_data pop (const event& ev)
    {
        const auto it = ev2data.find(ev);
        const auto evd = it->second;
        ev2data.erase(it);
        switch (evd.sc) {
        case SC_LATER:
            ev2t_later.erase(ev);
            _update_min_t_later(evd.t);
            break;
        case SC_SOONER:
            t2ev_sooner.erase(evd.t);
            break;
        case SC_NOW:
            ev_t_now.erase(_find_now(ev));
            n_ev_now--;
            break;
        case SC_NEVER:
            break;
        }
        return evd;
    }
    inline void erase (const event& ev)
    {
        // pop it without returning the data:
        pop(ev);
    }
    inline pair<event, event_data> pop_min_t ()
    {
        event ev;
        if (n_ev_now > 0)
        {
            // find minimal time by linear search
            // (this should be fastest here since ev_t_now is a very small set)
            timepoint min_t_now = INFINITY;
            event* ev_;
            auto it = ev_t_now.begin(), it_min = it;
            for (; it != ev_t_now.end(); it++)
            {
                auto t = it->second;
                if (t < min_t_now)
                {
                    it_min = it;
                    min_t_now = t;
                    ev_ = &(it->first);
                }
            }
            ev = *ev_;
            // erase it:
            ev_t_now.erase(it_min);
            n_ev_now--;
        }
        else
        {
            if (min_t_later < t2ev_sooner.begin()->first)
            {
                // copy all events from ev2t_later to t2ev_sooner:
                // TODO: is there a more efficient bulk insertion method?
                for (const auto& [ev2, t2] : ev2t_later) {
                    t2ev_sooner[t2] = ev2;
                    ev2data[ev2].sc = SC_SOONER;
                }
                // clear ev2t_later:
                ev2t_later.clear();
                min_t_later = INFINITY;
            }
            // pop first element from t2ev_sooner:
            const auto it = t2ev_sooner.begin();
            ev = it->second;
            t2ev_sooner.erase(it);
        }
        // finally pop corresponding data and return both:
        auto evd = ev2data.at(ev);
        ev2data.erase(ev);
        return pair<event, event_data>(ev, evd);
    }
    inline void insert (
            const event& ev,
            event_data evd      ///< [in] a copy (!) of the data to store
            )
    {
        assert(ev2data.count(ev) == 0);
        const auto t = evd.t;
        // depending on t, store in appropriate container:
        if (t > max_t)
        {
            evd.sc = SC_NEVER;
        }
        else if (t > current_t)
        {
            evd.sc = SC_LATER;
            ev2t_later[ev] = t;
            if (t < min_t_later) min_t_later = t;
        }
        else
        {
            evd.sc = SC_NOW;
            ev_t_now.push_back(pair<event, timepoint>(ev, t));
            n_ev_now++;
        }
        // store data:
        ev2data[ev] = evd;
    }
    inline void update (const event& ev, event_data& new_evd)
    {
        // (basically a combination of erase and insert)
        const auto old_evd = ev2data.at(ev);
        timepoint old_t = old_evd.t, new_t = new_evd.t;
        switch (old_evd.sc) {
        case SC_LATER:
            if (new_t > max_t)
            {
                // move to never:
                ev2t_later.erase(ev);
                _update_min_t_later(old_t);
                new_evd.sc = SC_NEVER;
            }
            else if (new_t > current_t)
            {
                // stay in later, just update t:
                new_evd.sc = SC_LATER;
                ev2t_later[ev] = new_t;
                if (new_t < min_t_later) min_t_later = new_t;
            }
            else
            {
                // move to now:
                ev2t_later.erase(ev);
                _update_min_t_later(old_evd.t);
                new_evd.sc = SC_NOW;
                ev_t_now.push_back(pair<event, timepoint>(ev, new_t));
                n_ev_now++;
            }
            break;
        case SC_SOONER:
            if (new_t > max_t)
            {
                // move to never:
                t2ev_sooner.erase(old_t);
                new_evd.sc = SC_NEVER;
            }
            else if (new_t > current_t)
            {
                // move to later:
                t2ev_sooner.erase(old_t);
                new_evd.sc = SC_LATER;
                ev2t_later[ev] = new_t;
                if (new_t < min_t_later) min_t_later = new_t;
            }
            else
            {
                // move to now:
                t2ev_sooner.erase(old_t);
                new_evd.sc = SC_NOW;
                ev_t_now.push_back(pair<event, timepoint>(ev, new_t));
                n_ev_now++;
            }
            break;
        case SC_NEVER:
            if (new_t > max_t)
            {
                // stay in never:
                new_evd.sc = SC_NEVER;
            }
            else if (new_t > current_t)
            {
                // move to later:
                new_evd.sc = SC_LATER;
                ev2t_later[ev] = new_t;
                if (new_t < min_t_later) min_t_later = new_t;
            }
            else
            {
                // move to now:
                new_evd.sc = SC_NOW;
                ev_t_now.push_back(pair<event, timepoint>(ev, new_t));
                n_ev_now++;
            }
            break;
        case SC_NOW:
            if (new_t > max_t)
            {
                // move to never:
                ev_t_now.erase(_find_now(ev));
                n_ev_now--;
                new_evd.sc = SC_NEVER;
            }
            else if (new_t > current_t)
            {
                // move to later:
                ev_t_now.erase(_find_now(ev));
                n_ev_now--;
                new_evd.sc = SC_LATER;
                ev2t_later[ev] = new_t;
                if (new_t < min_t_later) min_t_later = new_t;
            }
            else
            {
                // stay in now, update t:
                new_evd.sc = SC_NOW;
                _find_now(ev)->second = new_t;
            }
            break;
        }
        // store new data:
        ev2data[ev] = new_evd;
    }
};
