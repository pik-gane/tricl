/*
 * probability.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_PROBABILITY_H
#define INC_PROBABILITY_H

#include <assert.h>
#include <math.h>
#include <random>

#include "data_model.h"
#include "global_variables.h"

using namespace std;

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

inline probability probunit2probability (probunit pu, double left_tail, double right_tail)
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

inline rate effective_rate (rate attempt_rate, probunit success_probunit, double left_tail, double right_tail)
{
    rate r = attempt_rate * probunit2probability(success_probunit, left_tail, right_tail);
    assert (r >= 0);
    return r;
}

#endif



