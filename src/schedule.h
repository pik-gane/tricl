/** The schedule of the simulation: a sum tree over the rates of all scheduled events.
 *
 *  \file
 *
 *  The simulation uses the "direct" method: the waiting time until the next event is drawn from the total
 *  rate, and the event is then chosen with probability proportional to its rate. For this, the rates of all
 *  scheduled events (the effective rate of an individually scheduled event, the upper bound
 *  attempt rate * max. success probability of a summary event) are the leaves of a complete binary tree whose
 *  inner nodes hold the sums of their children, so that the total is always available exactly (up to the
 *  rounding of the current sums, without drift), an update costs O(log n), and choosing a leaf with
 *  probability proportional to its weight is a descent from the root.
 *
 *  Events with infinite rate ("immediate" events) are not in the tree but in a separate list (see
 *  immediate_events in global_variables.h), from which the next one is drawn uniformly at random.
 */

// make sure this file is only included once:
#ifndef INC_SCHEDULE_H
#define INC_SCHEDULE_H

#include <cmath>
#include <vector>

#include "data_model.h"

/** A complete binary tree of sums with events as leaves.
 */
class rate_tree
{
    std::vector<double> tree;   ///< tree[1] = total, leaves at indices cap ... 2 cap - 1, tree[i] = tree[2 i] + tree[2 i + 1]
    std::vector<event> slot2ev; ///< the event stored in each leaf slot
    std::vector<char> used;     ///< whether a leaf slot is in use
    std::vector<int> free_slots;  ///< the unused leaf slots
    int cap = 0;                ///< no. of leaf slots (a power of two, or 0)
    int n = 0;                  ///< no. of used leaf slots

    /** Recompute the sums on the path from leaf index i to the root. */
    inline void update_path (int i)
    {
        for (i /= 2; i >= 1; i /= 2) tree[i] = tree[2 * i] + tree[2 * i + 1];
    }

    /** Double the capacity, keeping all leaves. */
    void grow ()
    {
        int new_cap = (cap == 0) ? 64 : 2 * cap;
        std::vector<double> new_tree(2 * new_cap, 0.0);
        for (int s = 0; s < cap; s++) new_tree[new_cap + s] = tree[cap + s];
        for (int i = new_cap - 1; i >= 1; i--) new_tree[i] = new_tree[2 * i] + new_tree[2 * i + 1];
        tree.swap(new_tree);
        slot2ev.resize(new_cap);
        used.resize(new_cap, 0);
        for (int s = new_cap - 1; s >= cap; s--) free_slots.push_back(s);  // (lower slots are handed out first)
        cap = new_cap;
    }

public:

    /** \returns the total weight (sum of all rates in the tree). */
    inline double total () const { return (cap > 0) ? tree[1] : 0.0; }

    /** \returns the no. of events in the tree. */
    inline int size () const { return n; }

    /** \returns the weight of a leaf slot. */
    inline double weight (int slot) const { return tree[cap + slot]; }

    /** \returns the event in a leaf slot. */
    inline const event& at (int slot) const { return slot2ev[slot]; }

    /** \returns whether a leaf slot is in use. */
    inline bool is_used (int slot) const { return (slot >= 0) && (slot < cap) && used[slot]; }

    /** Insert an event with the given weight.
     *  \returns its leaf slot
     */
    inline int add (const event& ev, double w)
    {
        assert (w >= 0 && w < INFINITY);
        if (free_slots.empty()) grow();
        int s = free_slots.back();
        free_slots.pop_back();
        slot2ev[s] = ev;
        used[s] = 1;
        n++;
        tree[cap + s] = w;
        update_path(cap + s);
        return s;
    }

    /** Change the weight of a leaf slot. */
    inline void set (int slot, double w)
    {
        assert (is_used(slot) && (w >= 0) && (w < INFINITY));
        tree[cap + slot] = w;
        update_path(cap + slot);
    }

    /** Remove the event in a leaf slot. */
    inline void remove (int slot)
    {
        assert (is_used(slot));
        tree[cap + slot] = 0.0;
        update_path(cap + slot);
        used[slot] = 0;
        free_slots.push_back(slot);
        n--;
    }

    /** Choose the leaf in whose weight interval the number u falls, i.e. with probability proportional to its
     *  weight if u is uniformly distributed in [0, total).
     *  \returns the event in that leaf
     */
    inline const event& pick (double u) const
    {
        assert (n > 0 && total() > 0);
        int i = 1;
        while (i < cap) {
            i *= 2;
            if (u >= tree[i]) {
                u -= tree[i];
                i++;
            }
        }
        // guard against rounding: make sure we return a used leaf with positive weight
        if ((!used[i - cap]) || !(tree[i] > 0)) {
            int j = i;
            while ((j > cap) && ((!used[j - cap]) || !(tree[j] > 0))) j--;
            if ((!used[j - cap]) || !(tree[j] > 0)) {
                j = i;
                while ((j < 2 * cap - 1) && ((!used[j - cap]) || !(tree[j] > 0))) j++;
            }
            i = j;
        }
        return slot2ev[i - cap];
    }

    /** \returns the total recomputed from the leaves (for consistency checks). */
    double exact_total () const
    {
        double s = 0.0;
        for (int i = 0; i < cap; i++) if (used[i]) s += tree[cap + i];
        return s;
    }

    /** Remove all events. */
    void clear ()
    {
        tree.clear(); slot2ev.clear(); used.clear(); free_slots.clear();
        cap = n = 0;
    }
};

#endif
