/* (not yet used)
 *
 * Provide functionality for parameter estimation
 * via particle filtering with known (and hence common!) initial state
 *
 * Strategy
 * --------
 *
 * - replace current fast-access C arrays by efficient 2d C arrays
 * - replace current maps of numbers holding parameters by maps of C arrays
 * - allocate these via  ``new <type>[<length>]{}``
 * - in simulation mode, use length=1
 * - in pfilter mode:
 *   - use ev2evpfd (event_pfdata) instead of ev2evd and t2ev
 *   - use parse_event instead of pop_event and scheduling
 *   - this reads the next ec,e1,rat13,e3 from a csv file
 *   - initialize parameters according to some spec in config file
 * - combine resampling with mutation and/or gradient descent?
 *   - this needs derivatives of sigmoid function
 *
 *
 * Convention
 * ----------
 * p: particle
 */

struct event_pfdata {
    int n_angles = 0;                 ///< Current no. of angles influencing this event
    rate[] p2attempt_rate;            ///< Current attempt rate of this event by particle
    probunits[] p2success_probunits;  ///< Current success probunits of this event by particle
    rate[] p2effective_rate;          ///< Current attempt rate of this event by particle
};




