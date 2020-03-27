/*
 * probability.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

// for release mode, uncomment this:
#define NDEBUG
#include <assert.h>
//

#include <math.h>
#include <random>

#include "data_model.h"
#include "probability.h"

using namespace std;

// MATH:

// random generators:
random_device ran_dev;
mt19937 random_variable(ran_dev());
uniform_real_distribution<> uniform(0, 1);
exponential_distribution<> exponential(1);

double tail2scale(double tail)
{
    return 1 / (1 + 1 / tail) / pow(1 + log(1 + 1 / tail), tail + 1) / 2;
}

probability probunit2probability (probunit pu)
{
//    return 1.0 / (1.0 + exp(-l));
    double  left_tail = 1, // 1/sigma0 from paper. larger means heavier tail
            right_tail = 1,  // 1/sigma1
            scale = tail2scale(left_tail) + tail2scale(right_tail);
    return (1 / pow(1 + log(1 + exp(- pu / scale) / left_tail), left_tail)
            + 1
            - 1 / pow(1 + log(1 + exp(pu / scale) / right_tail), right_tail) ) / 2;
}

rate effective_rate (rate attempt_rate, probunit success_probunit)
{
    rate r = attempt_rate * probunit2probability(success_probunit);
    assert (r >= 0);
    return r;
}




