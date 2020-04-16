/** Some global constants.
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Mar 30, 2020
 *
 * @file
 */

#include "data_model.h"

namespace tricl {

    unordered_map<event_class, string> ec2label = {
        { EC_EST, "establish that" },
        { EC_TERM, "terminate that" },
        { EC_ACT, "let" }
    };

}

