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
extern mt19937 random_variable;                ///< Our pseudo-random number generator for the dynamics (scheduling of events, summary event draws)
extern mt19937 initial_state_random_variable;  ///< A separate generator for the random initial state (block and geometric models), so that the initial state depends only on the config and the seed
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
        double right_tail, ///< [in] the tail index of the right (upper) tail, >= 0
        double scale       ///< [in] the precomputed scale parameter, = tail2scale(left_tail) + tail2scale(right_tail)
        )
{
    if ((left_tail == 0) && (right_tail == 0))
    {
        return 1 / (1 + exp(- pu));
    }
    else
    {
        return (      pow(1 + log(1 + left_tail  * exp(- pu / scale)),
                          - 1 / left_tail)
                + 1 - pow(1 + log(1 + right_tail * exp(  pu / scale)),
                          - 1 / right_tail)
               ) / 2;
    }
    // TODO: what if only one tail index == 0 ?
}

/** Convert probability units to probability, computing the scale parameter on the fly
 *  (for non-performance-critical uses; the simulation uses the precomputed scale from evtid2scale).
 */
inline probability probunits2probability (probunits pu, double left_tail, double right_tail)
{
    return probunits2probability(pu, left_tail, right_tail, tail2scale(left_tail) + tail2scale(right_tail));
}

/** Derivative of \ref probunits2probability() w.r.t. the probability units.
 *
 *  \returns the derivative, >= 0 (0 for infinite probability units)
 */
inline double probunits2probability_derivative (
        probunits pu,      ///< [in] the probability units, -inf...inf
        double left_tail,  ///< [in] the tail index of the left (lower) tail, >= 0
        double right_tail, ///< [in] the tail index of the right (upper) tail, >= 0
        double scale       ///< [in] the precomputed scale parameter, = tail2scale(left_tail) + tail2scale(right_tail)
        )
{
    if (!std::isfinite(pu)) return 0.0;
    if ((left_tail == 0) && (right_tail == 0))
    {
        double p = 1 / (1 + exp(- pu));
        return p * (1 - p);
    }
    else
    {
        // d/dpu of pow(1 + log(1 + L * exp(-pu/s)), -1/L) = v^(-1/L-1) * exp(-pu/s) / (s * u) with u = 1 + L exp(-pu/s), v = 1 + log(u);
        // (this mirrors the formula in probunits2probability, including its behaviour for a single zero tail index)
        double left = 0.0, right = 0.0;
        if (left_tail > 0) {
            double e = exp(- pu / scale);
            if (std::isfinite(e)) {
                double u = 1 + left_tail * e, v = 1 + log(u);
                left = pow(v, - 1 / left_tail - 1) * e / (scale * u);
            }
        }
        if (right_tail > 0) {
            double e = exp(pu / scale);
            if (std::isfinite(e)) {
                double u = 1 + right_tail * e, v = 1 + log(u);
                right = pow(v, - 1 / right_tail - 1) * e / (scale * u);
            }
        }
        return (left + right) / 2;
    }
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
        double right_tail,      ///< [in] the right tail index of the sigmoid function to be used, >= 0
        double scale            ///< [in] the precomputed scale parameter of the sigmoid function
        )
{
    assert (attempt_rate >= 0);
    if (attempt_rate == 0) return 0;
    if (success_pus == -INFINITY) return 0;  // impossible events never happen, even if attempted at an infinite rate
    if (attempt_rate == INFINITY) return INFINITY;
    rate r = attempt_rate * probunits2probability(success_pus, left_tail, right_tail, scale);
    assert (r >= 0);
    return r;
}

/** Compute the effective rate, computing the sigmoid's scale parameter on the fly (for non-performance-critical uses).
 */
inline rate effective_rate (rate attempt_rate, probunits success_pus, double left_tail, double right_tail)
{
    return effective_rate(attempt_rate, success_pus, left_tail, right_tail, tail2scale(left_tail) + tail2scale(right_tail));
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
 *
 *  Since the total is maintained incrementally in floating point arithmetic,
 *  it may drift by a tiny amount; hence only clearly negative values are treated as an inconsistency
 *  (and only if do_assert is set: when rescheduling a summary event, the total may legitimately become
 *  negative for a moment, since the summary event's rate also covers pairs whose share has been subtracted).
 *  The total is also recomputed exactly every now and then in step().
 */
inline void subtract_effective_rate (rate er, bool do_assert)
{
    if (er < INFINITY)
    {
        total_finite_effective_rate -= er;
        if (debug) cout << "             finite er - " << er << " = " << total_finite_effective_rate << endl;
        if (do_assert) assert (total_finite_effective_rate > -1e-6 * max(1.0, er));
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
