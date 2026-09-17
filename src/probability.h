// make sure this file is only included once:
#ifndef INC_PROBABILITY_H
#define INC_PROBABILITY_H

/** Inline functions dealing with probabilities
 *
 *  \file
 */

#include <cmath>
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

const double ln2 = std::log(2.0);  ///< ln 2, used in tail2scale()

void init_randomness ();

/** Softplus function ln(1 + e^x), evaluated without overflow.
 *
 *  Auxiliary function for \ref tail_term().
 */
inline double softplus (
        double x  ///< [in] any real number, may be +-inf
        )
{
    return (x > 0) ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

/** Expit (inverse logit) function 1 / (1 + e^-x).
 */
inline double expit (
        double x  ///< [in] any real number, may be +-inf
        )
{
    return 1 / (1 + std::exp(-x));
}

/** Tail term T_tail(x) = (1 + tail * ln(1 + e^x))^(-1/tail) of the sigmoidal function.
 *
 *  For tail == 0 this is its limit 1 / (1 + e^x) (the mirrored expit function),
 *  for tail > 0 the exponential exp(-y), y = ln(1 + e^x), of that expit function
 *  is replaced by the power-law type function (1 + tail * y)^(-1/tail),
 *  so that T_tail(x) ~ (tail * x)^(-1/tail) as x -> +inf,
 *  while 1 - T_tail(x) ~ e^x as x -> -inf for all tail indices.
 *  T_tail(x) decreases from 1 to 0 and depends continuously on tail, including at tail == 0.
 *
 *  Auxiliary function for \ref probunits2probability().
 *
 *  \returns the tail term, 0...1
 */
inline double tail_term (
        double tail,  ///< [in] tail index, >= 0
        double x      ///< [in] argument, -inf...inf
        )
{
    double y = softplus(x);
    return (tail == 0)
            ? std::exp(- y)
            : std::exp(- std::log1p(tail * y) / tail);
}

/** Complement 1 - T_tail(x) of the tail term, evaluated without cancellation for x -> -inf.
 *
 *  Auxiliary function for \ref probunits2probability().
 *
 *  \returns 1 - tail_term(tail, x), 0...1
 */
inline double tail_term_complement (
        double tail,  ///< [in] tail index, >= 0
        double x      ///< [in] argument, -inf...inf
        )
{
    double y = softplus(x);
    return (tail == 0)
            ? - std::expm1(- y)
            : - std::expm1(- std::log1p(tail * y) / tail);
}

/** Negative derivative -dT_tail(x)/dx = (1 + tail * ln(1 + e^x))^(-1/tail - 1) * expit(x) of the tail term.
 *
 *  Auxiliary function for \ref probunits2probability_derivative().
 *
 *  \returns the negative derivative, >= 0
 */
inline double tail_term_slope (
        double tail,  ///< [in] tail index, >= 0
        double x      ///< [in] argument, -inf...inf
        )
{
    double y = softplus(x);
    double base = (tail == 0)
            ? std::exp(- y)
            : std::exp(- (1 + 1 / tail) * std::log1p(tail * y));
    return base * expit(x);
}

/** Compute the scale factor k(tail) = (1 + tail * ln 2)^(-1 - 1/tail) for a tail index.
 *
 *  k(0) = 1/2 is the limit for tail -> 0.
 *  The sum of the two scale factors of an event type normalises the probability units
 *  so that the sigmoidal function has slope 1/4 at zero probability units for all tail indices,
 *  the same slope as the expit function.
 *
 *  Auxiliary function for \ref probunits2probability().
 *
 *  \returns the scale factor, 0 < k <= 1/2
 */
inline double tail2scale (
        double tail  ///< [in] tail index to compute the scale factor for, >= 0
        )
{
    return (tail == 0)
            ? 0.5
            : std::exp(- (1 + 1 / tail) * std::log1p(tail * ln2));
}

/** Convert probability units to probability.
 *
 *  This is a smooth sigmoidal function
 *  that also depends continuously on two tail indices as parameters:
 *
 *      f(pu) = T_left(-v) / 2 + 1/2 - T_right(v) / 2,   v = pu / (k(left) + k(right)),
 *
 *  with the tail term T (see \ref tail_term()) and the scale factors k (see \ref tail2scale()).
 *  If both tail indices are zero, this is exactly the expit (= inverse logit) function 1 / (1 + e^-pu).
 *  If a tail index is positive, the corresponding tail
 *  converges with a power-law decay to its limit 0 (left tail) or 1 (right tail),
 *  where the power-law exponent is 1 / tail index;
 *  if it is zero, the corresponding tail converges exponentially.
 *  Any combination of tail indices is allowed, including exactly one zero tail index.
 *  For all tail indices the function is strictly increasing and has slope 1/4 at pu == 0.
 *
 *  \returns the probability, 0...1
 */
inline probability probunits2probability (
        probunits pu,      ///< [in] the probability units, -inf...inf
        double left_tail,  ///< [in] the tail index of the left (lower) tail, >= 0
        double right_tail, ///< [in] the tail index of the right (upper) tail, >= 0
        double scale       ///< [in] the precomputed scale, = tail2scale(left_tail) + tail2scale(right_tail)
        )
{
    if ((left_tail == 0) && (right_tail == 0))
    {
        // (this is what the general formula gives for these tail indices, computed more cheaply)
        return expit(pu);
    }
    double v = pu / scale;
    // both summands vanish for pu -> -inf, so small probabilities are computed without cancellation:
    return (tail_term(left_tail, - v) + tail_term_complement(right_tail, v)) / 2;
}

/** Convert probability units to probability, computing the scale on the fly
 *  (for non-performance-critical uses; the simulation uses the precomputed scale from evtid2scale).
 */
inline probability probunits2probability (probunits pu, double left_tail, double right_tail)
{
    return probunits2probability(pu, left_tail, right_tail, tail2scale(left_tail) + tail2scale(right_tail));
}

/** Derivative of \ref probunits2probability() w.r.t. the probability units.
 *
 *  \returns the derivative, >= 0 (0 for infinite probability units); 1/4 at pu == 0 for all tail indices
 */
inline double probunits2probability_derivative (
        probunits pu,      ///< [in] the probability units, -inf...inf
        double left_tail,  ///< [in] the tail index of the left (lower) tail, >= 0
        double right_tail, ///< [in] the tail index of the right (upper) tail, >= 0
        double scale       ///< [in] the precomputed scale, = tail2scale(left_tail) + tail2scale(right_tail)
        )
{
    if (!std::isfinite(pu)) return 0.0;
    if ((left_tail == 0) && (right_tail == 0))
    {
        double p = expit(pu);
        return p * (1 - p);
    }
    double v = pu / scale;
    return (tail_term_slope(left_tail, - v) + tail_term_slope(right_tail, v)) / (2 * scale);
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
