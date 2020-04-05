/*
 * debugging.cpp
 *
 *  Created on: Apr 5, 2020
 *      Author: heitzig
 */

#include "assert.h"
#include <iostream>
#include <string>

#include "global_variables.h"
#include "io.h"

void verify_data_consistency() {
    // ev2data:
    for (auto& [ev, evd] : ev2data) {
        assert (evd.n_angles >= 0);
        assert (evd.attempt_rate > 0);
        assert (evd.success_probunits > -INFINITY);
        if (!(evd.t > -INFINITY)) dump_data();
        assert (evd.t > -INFINITY);
        assert (t2be.count(evd.t) == 1);
    }
    // t2be:
    for (auto& [t, ev] : t2be) {
        assert (t > -INFINITY);
        assert (ev2data.count(ev) == 1);
    }
}
