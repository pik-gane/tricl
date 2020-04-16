/*
 * probability.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_PROBABILITY_H
#define INC_PROBABILITY_H

#include <random>

#include "global_variables.h"

using std::mt19937;
using std::random_device;
using std::uniform_real_distribution;
using std::exponential_distribution;

// random generators:
extern mt19937 random_variable;
extern uniform_real_distribution<> uniform;
extern exponential_distribution<> exponential;

const double scale0 =  1 / 2 / exp(1);

void init_randomness ();

inline double tail2scale(double tail)
{
    return (tail == 0) ? scale0 : 1 / (1 + tail) / pow(1 + log(1 + tail), 1 + 1 / tail) / 2;
}

inline probability probunits2probability (probunits pu, double left_tail, double right_tail)
{
    if ((left_tail == 0) && (right_tail == 0)) {
        return 1 / (1 + exp(-pu));
    } else {
        double scale = tail2scale(left_tail) + tail2scale(right_tail);
        return (1 / pow(1 + log(1 + left_tail * exp(- pu / scale)), 1 / left_tail)
                + 1
                - 1 / pow(1 + log(1 + right_tail * exp(pu / scale)), 1 / right_tail) ) / 2;
    }
}

inline rate effective_rate (rate attempt_rate, probunits success_pus, double left_tail, double right_tail)
{
    assert (attempt_rate >= 0);
    if (attempt_rate == 0) return 0;
    if (attempt_rate == INFINITY) return INFINITY;
    rate r = attempt_rate * probunits2probability(success_pus, left_tail, right_tail);
    assert (r >= 0);
    return r;
}

#endif



