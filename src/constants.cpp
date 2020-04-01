/** Some global constants.
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Mar 30, 2020
 *
 * @file
 */

#include "data_model.h"

unordered_map<event_class, string> ec2label = {
    { EC_EST, "establish" },
    { EC_TERM, "terminate" },
    { EC_ACT, "act" }
};



