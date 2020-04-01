/*
 * probability.h
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#ifndef INC_PROBABILITY_H
#define INC_PROBABILITY_H

#include <random>

#include "data_model.h"

using namespace std;

// MATH:

// random generators:
extern random_device ran_dev;
extern mt19937 random_variable;
extern uniform_real_distribution<> uniform;
extern exponential_distribution<> exponential;

double tail2scale(double tail);

probability probunit2probability (probunit pu);

rate effective_rate (rate attempt_rate, probunit success_probunit);

#endif



