/** Unit tests of the sigmoidal function probunits2probability() and its derivative (see src/probability.h).
 *
 *  \file
 *
 *  Reference values are from the handover note that introduced the q-exponential form of the function
 *  (which was normalised to slope 1/4 at zero; the function is now normalised to slope 1, i.e. its argument is
 *  multiplied by 4, so the reference values are checked at a quarter of the handover's probability units).
 *  Prints the failed checks and returns a non-zero exit code if any check fails.
 */

#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include "probability.h"

static int n_failed = 0, n_checks = 0;

static void check (bool ok, const std::string& what)
{
    n_checks++;
    if (!ok) {
        n_failed++;
        std::printf("FAILED: %s\n", what.c_str());
    }
}

static bool near_abs (double a, double b, double tol) { return std::fabs(a - b) <= tol; }
static bool near_rel (double a, double b, double rel) { return std::fabs(a - b) <= rel * std::fabs(b); }

static double slope_fd (double pu, double l, double r)
{
    double h = 1e-6 * std::max(1.0, std::fabs(pu));
    return (probunits2probability(pu + h, l, r) - probunits2probability(pu - h, l, r)) / (2 * h);
}

int main ()
{
    char buf[200];

    // reference values (at a quarter of the handover's probability units, see above):
    check(near_abs(probunits2probability(0.25, 0, 0), 0.731058578630, 1e-11), "(0,0) pu=1/4 -> expit(1)");
    check(near_rel(probunits2probability(-12.5, 0, 1), 2.6e-26, 0.05), "(0,1) pu=-50/4 -> 2.6e-26");
    check(near_abs(probunits2probability(12.5, 0, 1), 0.991653, 1e-6), "(0,1) pu=+50/4 -> 0.991653");
    check(near_abs(probunits2probability(0.75, 1e-9, 1), 0.876248035533, 1e-10), "(1e-9,1) pu=3/4");
    check(near_abs(probunits2probability(0.75, 0, 1), 0.876248035498, 1e-10), "(0,1) pu=3/4");
    check(near_abs(probunits2probability(0.75, 1e-9, 1) - probunits2probability(0.75, 0, 1), 0.0, 1e-8), "continuity at zero tail index");

    // for the default tail indices (1,1), the function equals the one of earlier versions of tricl,
    // f(pu) = ((1 + ln(1 + e^-pu/s))^-1 + 1 - (1 + ln(1 + e^pu/s))^-1) / 2 with s = 1 / (2 (1 + ln 2)^2):
    {
        double s = 1 / (2 * std::pow(1 + std::log(2.0), 2));
        for (double pu = -20; pu <= 20; pu += 0.125) {
            double old = (std::pow(1 + std::log(1 + std::exp(- pu / s)), -1.0) + 1 - std::pow(1 + std::log(1 + std::exp(pu / s)), -1.0)) / 2;
            std::snprintf(buf, sizeof buf, "earlier formula for tails (1,1) at pu=%g: %.17g vs %.17g", pu, probunits2probability(pu, 1, 1), old);
            check(near_abs(probunits2probability(pu, 1, 1), old, 1e-14), buf);
        }
    }

    // slope 1 at zero probability units, by finite differences and analytically:
    for (auto tails : std::vector<std::pair<double, double>> { {0, 0}, {1, 1}, {2, 0.5}, {0, 1}, {1, 0}, {5, 0.2}, {10, 10} }) {
        double l = tails.first, r = tails.second;
        std::snprintf(buf, sizeof buf, "slope at 0 for tails (%g,%g)", l, r);
        check(near_abs(slope_fd(0, l, r), 1.0, 1e-6), buf);
        std::snprintf(buf, sizeof buf, "analytic derivative at 0 for tails (%g,%g)", l, r);
        check(near_abs(probunits2probability_derivative(0, l, r, tail2scale(l) + tail2scale(r)), 1.0, 1e-12), buf);
        // limits:
        std::snprintf(buf, sizeof buf, "limits for tails (%g,%g)", l, r);
        check(probunits2probability(INFINITY, l, r) == 1.0 && probunits2probability(-INFINITY, l, r) == 0.0, buf);
        check(probunits2probability_derivative(INFINITY, l, r, tail2scale(l) + tail2scale(r)) == 0.0, buf);
    }

    // power-law lower tail with exponent -1/tail for tail index 1:
    {
        double e = std::log(probunits2probability(-1e4, 1, 1) / probunits2probability(-1e3, 1, 1)) / std::log(10.0);
        std::snprintf(buf, sizeof buf, "lower tail exponent for tails (1,1): %.4f (expected -1)", e);
        check(near_abs(e, -1.0, 0.01), buf);
    }

    // exactly the expit function of 4 pu for tail indices (0,0), also when computed by the general formula:
    for (double pu = -30; pu <= 30; pu += 0.5) {
        std::snprintf(buf, sizeof buf, "expit at pu=%g", pu);
        check(near_rel(probunits2probability(pu, 0, 0), 1 / (1 + std::exp(-4 * pu)), 1e-14), buf);
        double v = pu / (tail2scale(0) + tail2scale(0));
        check(near_rel((tail_term(0, -v) + tail_term_complement(0, v)) / 2, 1 / (1 + std::exp(-4 * pu)), 1e-13), buf);
    }

    // monotonicity and derivative on a grid of tail indices:
    std::vector<double> tail_grid { 0, 0.1, 0.5, 1, 2, 5, 10 };
    for (double l : tail_grid) {
        for (double r : tail_grid) {
            double scale = tail2scale(l) + tail2scale(r), last = -1;
            for (double pu = -100; pu <= 100; pu += 0.25) {
                double p = probunits2probability(pu, l, r, scale);
                std::snprintf(buf, sizeof buf, "range/monotonicity for tails (%g,%g) at pu=%g: %g after %g", l, r, pu, p, last);
                check(std::isfinite(p) && p >= 0 && p <= 1 && p >= last, buf);
                if ((pu > -5) && (pu <= 5)) check(p > last, buf);  // (strictly increasing where not saturated in double precision)
                last = p;
                double d = probunits2probability_derivative(pu, l, r, scale);
                std::snprintf(buf, sizeof buf, "derivative for tails (%g,%g) at pu=%g: %g vs finite difference %g", l, r, pu, d, slope_fd(pu, l, r));
                check(d >= 0 && near_abs(d, slope_fd(pu, l, r), 1e-7 + 1e-6 * d), buf);
            }
        }
    }

    // the sigmoid is 1/2 at zero probability units for equal tail indices:
    for (double t : tail_grid) {
        std::snprintf(buf, sizeof buf, "f(0) = 1/2 for tails (%g,%g)", t, t);
        check(near_abs(probunits2probability(0, t, t), 0.5, 1e-14), buf);
    }

    std::printf("%d of %d checks passed\n", n_checks - n_failed, n_checks);
    return n_failed ? 1 : 0;
}
