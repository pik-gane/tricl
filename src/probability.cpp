/** Additional functions dealing with probabilities.
 *
 *  \file
 */

#include "global_variables.h"
#include "probability.h"

// random generators:
random_device ran_dev;
mt19937 random_variable;
uniform_real_distribution<> uniform(0, 1);
exponential_distribution<> exponential(1);

/** Initialize the pseudo-random number generator using specified seed.
 *
 *  Uses a pseudo-random seed if seed == 0.
 */
void init_randomness ()
{
    auto theseed = (seed == 0) ? ran_dev() : seed;
    if (!quiet) cout << " using random seed " << theseed << endl;
    random_variable = mt19937(theseed);
}




