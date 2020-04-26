// make sure this file is only included once:
#ifndef INC_PROBABILITY_H
#define INC_PROBABILITY_H

/** Inline functions dealing with probabilities
 *
 *  \file
 */

#include <random>

#include "global_variables.h"

using std::mt19937;
using std::random_device;
using std::uniform_real_distribution;
using std::exponential_distribution;

// random generators:
extern mt19937 random_variable;                 ///< Our pseudo-random number generator
extern uniform_real_distribution<> uniform;     ///< uniform(random_variable) produces uniformly distributed numbers 0...1
extern exponential_distribution<> exponential;  ///< exponential(random_variable) produces exponentially distributed numbers with mean 1

const double scale0 =  1 / 2 / exp(1);  ///< Precomputed scale parameter for tail index 0

void init_randomness ();

/** Compute the scale parameter for a tail index.
 *
 *  Auxiliary function for \ref probunits2probability().
 *
 *  \returns the scale parameter, > 0
 */
inline double tail2scale (
        double tail  ///< [in] tail index to compute the scale parameter for, >= 0
        )
{
    return (tail == 0)
            ? scale0
            : 1 / (1 + tail)
                / pow(1 + log(1 + tail),
                      1 + 1 / tail)
                / 2;
}

/** Convert probability units to probability.
 *
 *  This is a smooth sigmoidal function
 *  that also depends smoothly on two tail indices as parameters.
 *  If both tail indices are zero, this is simply the expit (= inverse logit) function.
 *  If a tail index is positive, the corresponding tail
 *  converges with a power-law decay to its limit 0 (left tail) or 1 (right tail),
 *  where the power-law exponent is 1 / tail index.
 *
 *  \returns the probability, 0...1
 */
inline probability probunits2probability (
        probunits pu,      ///< [in] the probability units, -inf...inf
        double left_tail,  ///< [in] the tail index of the left (lower) tail, >= 0
        double right_tail  ///< [in] the tail index of the right (upper) tail, >= 0
        )
{
    if ((left_tail == 0) && (right_tail == 0))
    {
        return 1 / (1 + exp(- pu));
    }
    else
    {
        double scale = tail2scale(left_tail) + tail2scale(right_tail);

        return (      pow(1 + log(1 + left_tail  * exp(- pu / scale)),
                          - 1 / left_tail)
                + 1 - pow(1 + log(1 + right_tail * exp(  pu / scale)),
                          - 1 / right_tail)
               ) / 2;
    }
    // TODO: what if only one tail index == 0 ?
}

/** Compute the current effective rate at which an event occurs
 *  from its current attempt rate and success probability units.
 *
 *  An effective rate of inf implies that the event occurs immediately.
 *
 *  \returns the effective rate, 0...inf
 */
inline rate effective_rate (
        rate attempt_rate,      ///< [in] the event's current total attempt rate, 0...inf
        probunits success_pus,  ///< [in] the event's current total success probability units, -inf...inf
        double left_tail,       ///< [in] the left tail index of the sigmoid function to be used, >= 0
        double right_tail       ///< [in] the right tail index of the sigmoid function to be used, >= 0
        )
{
    assert (attempt_rate >= 0);
    if (attempt_rate == 0) return 0;
    if (attempt_rate == INFINITY) return INFINITY;  // even if pus == -inf !
    rate r = attempt_rate * probunits2probability(success_pus, left_tail, right_tail);
    assert (r >= 0);
    return r;
}

/** Add effective rate to total, taking care of infinite values.
 */
inline void add_effective_rate (rate er)
{
    if (er < INFINITY) {
        total_finite_effective_rate += er;
        if (debug) cout << "             finite er + " << er << " = " << total_finite_effective_rate << endl;
    }
    else {
        n_infinite_effective_rates++;
        if (debug) cout << "             infinite ers + 1 = " << n_infinite_effective_rates << endl;
    }
}

/** Subtract effective rate to total, taking care of infinite values.
 */
inline void subtract_effective_rate (rate er, bool do_assert)
{
    if (er < INFINITY)
    {
        total_finite_effective_rate -= er;
        if (debug) cout << "             finite er - " << er << " = " << total_finite_effective_rate << endl;
        if (do_assert) assert (total_finite_effective_rate >= 0);
    }
    else
    {
        n_infinite_effective_rates--;
        if (debug) cout << "             infinite ers - 1 = " << n_infinite_effective_rates << endl;
        if (do_assert) assert (n_infinite_effective_rates >= 0);
    }
}
inline void subtract_effective_rate (rate er)
{
    subtract_effective_rate(er, true);
}

/** \returns actual total effective rate, taking care of infinite values.
 */
inline rate total_effective_rate ()
{
    return (n_infinite_effective_rates > 0) ? INFINITY : total_finite_effective_rate;
}


#endif
