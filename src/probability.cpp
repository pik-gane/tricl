/*
 * probability.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include "global_variables.h"
#include "probability.h"

// random generators:
random_device ran_dev;
mt19937 random_variable;
uniform_real_distribution<> uniform(0, 1);
exponential_distribution<> exponential(1);

void init_randomness ()
{
    auto theseed = (seed == 0) ? ran_dev() : seed;
    if (!quiet) cout << " using random seed " << theseed << endl;
    random_variable = mt19937(theseed);
}




