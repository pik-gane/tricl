/** Some global constants.
 *
 *  \file
 */

#include "data_model.h"


namespace tricl {

    /** Labels for event classes.
     *
     *  These labels are prepended to the event label
     *  (which is a simple sentence like "A loves B").
     */
    unordered_map<event_class, string> ec2label = {
        { EC_EST, "establish that" },    ///< establishment of a relationship
        { EC_TERM, "terminate that" },   ///< termination of a relationship
        { EC_ACT, "let it occur that" }  ///< occurrence of an action
    };

} // end of namespace tricl


/** Event data for initial events that lead to logl = 0 for these events
 */
event_data sure_evd = {
        .n_angles = 0,
        .attempt_rate = 0,
        .n_inf_attempt = 1,        // infinite attempt rate
        .success_probunits = 0,
        .n_pos_inf_probunits = 1,  // certain success
        .n_neg_inf_probunits = 0,
        .effective_rate = INFINITY,
        .t = -INFINITY
};

