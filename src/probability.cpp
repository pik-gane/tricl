/** Additional functions dealing with probabilities.
 *
 *  \file
 */

#include "global_variables.h"
#include "probability.h"

// random generators:
random_device ran_dev;
mt19937 random_variable;
mt19937 initial_state_random_variable;
uniform_real_distribution<> uniform(0, 1);
exponential_distribution<> exponential(1);

/** Initialize the pseudo-random number generators using the specified seed.
 *
 *  Uses a pseudo-random seed if seed == 0.
 *  The generator for the random initial state is seeded differently from the one for the dynamics,
 *  so that the initial state depends only on the config file and the seed (and not on whether the
 *  dynamics are simulated or replayed).
 */
void init_randomness ()
{
    if (seed == 0) seed = ran_dev();  // store the actually used seed for reporting
    if (!quiet) cout << " using random seed " << seed << endl;
    random_variable = mt19937(seed);
    std::seed_seq initial_state_seed { seed, 1u };
    initial_state_random_variable = mt19937(initial_state_seed);
}




