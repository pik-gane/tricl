// make sure this file is only included once:
#ifndef INC_ANGLE_H
#define INC_ANGLE_H

/** Performance-critical inline functions for handling of angles.
 *
 *  \file
 *
 *  See \ref data_model.h for how a \ref tricl::angle relates to other tricl datatypes.
 */

#include "global_variables.h"
#include "probability.h"
#include "debugging.h"
#include "event.h"
#include "io.h"

/** Perform all necessary changes in state and event data
 *  due to the addition or deletion of an angle.
 *
 *  Iterate through all events that might be influenced by this angle,
 *  update their attempt rates and success probability units,
 *  and (re)schedule them based on the new rates.
 *
 *  This is one of the performance bottleneck functions
 *  since it is called many times by \ref update_adjacent_events().
 */
inline void add_or_delete_angle (
        event_class ec_angle,               ///< [in] class of the event that caused the change: EC_EST adds an angle, EC_TERM deletes one
        entity e1,                          ///< [in] source entity
        entity_type et1,                    ///< [in] its type
        relationship_or_action_type rat12,  ///< [in] source-to-middle relationship or action type
        entity e2,                          ///< [in] middle entity
        entity_type et2,                    ///< [in] its type
        relationship_or_action_type rat23,  ///< [in] middle-to-target relationship or action type
        entity e3,                          ///< [in] target entity
        entity_type et3                     ///< [in] its type
        )
{
    if (debug) cout << "    " << ec2label[ec_angle] << " \"" << e2label[e1] << " " << rat2label[rat12] << " "
            << e2label[e2] << " " << rat2label[rat23] << " " << e2label[e3] << "\"" << endl;

    // update total no. of angles (if non-id.):
    n_angles += ((e1 == e2) || (e2 == e3) || (e3 == e1)) ? 0 : (ec_angle == EC_EST) ? 1 : -1;

    // iterate through all possible source-target relationship or action types:
    for (auto& rat13 : ets2relations[{ et1, et3 }]) {
        bool link13_exists = (e2outs[e1].count({ .rat_out = rat13, .e_target = e3 }) > 0);

        // construct the type of the corresponding event whose data might need an update:
        event_class ec13 = link13_exists ? EC_TERM : EC_EST;
        event_type evt = { .ec=ec13, et1, rat13, et3 };
        if (debug) cout << "     possibly updating event: " << ec2label[ec13] <<  " \"" << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << "\"" << endl;

        // only continue if the event type can happen at all:
        if (possible_evts.count(evt) > 0) {
            event ev = { .ec=ec13, e1, rat13, e3 };
            if (debug) cout << "      event type " << evt << " has a base success prob. of "
                    << probunits2probability(evt2base_probunits.at(evt), evt2left_tail.at(evt), evt2right_tail.at(evt)) << endl;

            // get influence of angle on event:
            influence_type inflt = { .evt = evt, .at = { rat12, et2, rat23 } };
            auto dar = _inflt2attempt_rate[INFLT(inflt)];
            auto dspu = _inflt2delta_probunits[INFLT(inflt)];

            // only continue if influence is nonzero:
            if (COUNT_ALL_ANGLES || (dar != 0.0) || (dspu != 0.0)) {
                if (debug) cout << "       angle may influence attempt or success" << endl;
                if (ec_angle == EC_EST) { // angle is added:
                    if (ev2data.count(ev) == 0) {
                        if (debug) cout << "        event will be scheduled newly" << endl;
                        auto evd_ = &ev2data[ev];  // generates a new event_data object
                        evd_->n_angles = 1;
                        evd_->attempt_rate = evt2base_attempt_rate.at(evt) + dar;
                        evd_->success_probunits = evt2base_probunits.at(evt) + dspu;
                        schedule_event(ev, evd_, evt2left_tail.at(evt), evt2right_tail.at(evt));
                    } else {
                        if (debug) cout << "        event will be rescheduled" << endl;
                        auto evd_ = &(ev2data.at(ev));
                        evd_->n_angles += 1;
                        evd_->attempt_rate += dar;
                        evd_->success_probunits += dspu;
                        reschedule_event(ev, evd_, evt2left_tail.at(evt), evt2right_tail.at(evt));
                    }
                } else { // angle is removed:
                    auto evd_ = &(ev2data.at(ev));
                    assert (evd_->n_angles > 0);  // since angle must have been added earlier to be removed now
                    evd_->n_angles -= 1;
                    evd_->attempt_rate = max(0.0, evd_->attempt_rate - dar);
                    evd_->success_probunits -= dspu;
                    if ((ec13 != EC_TERM) && (evd_->n_angles == 0)) { // only spontaneous non-termination event is left:
                        if (debug) cout << "        last angle was removed, so event will be removed" << endl;
                        // remove specific event:
                        remove_event(ev, evd_);  // event must have been scheduled earlier when angle was added
                    } else {
                        if (debug) cout << "        event will be rescheduled" << endl;
                        reschedule_event(ev, evd_, evt2left_tail.at(evt), evt2right_tail.at(evt));  // event must have been scheduled earlier when angle was added
                    }
                }
            } else {
                if (debug) cout << "       but angle may not influence attempt or success" << endl;
            }
        }
    }
}

/** Compare each \ref outleg of e1 with each \ref inleg of e3 to find each \ref angle from e1 to e3.
 *
 *  This is one of the performance bottleneck functions
 *  since it is called by \ref add_or_delete_angle().
 *  It uses a large share of the model's CPU time.
 *
 *  (code was adapted from adapted from set_intersection template)
 *
 *  \returns a vector of found angles
 */
inline angle_vec leg_intersection (
        const entity e1,         ///< [in] source entity
        const outleg_set& out1,  ///< [in] set of outlegs of source entity
        const inleg_set& in3,    ///< [in] set of inlegs of target entity
        const entity e3          ///< [in] target entity
        )
{
  // allocate mem for result:
  angle_vec result(max(out1.size(), in3.size()) * n_rats * n_rats);
  // work with pointers (iterators):
  auto out1_it = out1.begin();
  auto in3_it = in3.begin();
  auto result_it = result.begin();
  /** Algorithm:
   *  ----------
   *  The two sequences are sorted by e2 (since std::set is an ordered datatype and operator< for legs was implemented accordingly).
   *  Pseudocode:
   *
   *      put blockstart = out1.end().
   *      repeat:
   *        if out1.e2 < in3.e2, advance out1.
   *        else if out1.e2 > in3.e2:
   *          advance in3.
   *          if blockstart != out1.end():
   *            if in3.e2 == previous in3.e2, rewind out2 to blockstart
   *            else put blockstart = out1.end().
   *        else out1.e2 == in3.e2:
   *          if blockstart == out1.end(), remember out1 position as blockstart
   *          store found angle
   *          advance out1
   */
  auto out1end = out1.end(), blockstart = out1end;
  auto in3end = in3.end();
  entity last_e2 = out1_it->e_target;
  while (in3_it != in3end)
  {
      if (out1_it == out1end) {
          if (blockstart == out1end) {
              break;
          } else {
              ++in3_it;
              if (in3_it == in3end) {
                  break;
              }
              if (in3_it->e_source == last_e2) {
                  out1_it = blockstart;
              } else {
                  break;
              }
          }
      }
      if (debug) cout << "         checking: " << rat2label[out1_it->rat_out] << " " << e2label[out1_it->e_target] << ", "
              << e2label[in3_it->e_source] << " " << rat2label[in3_it->rat_in] << endl;
      if (out1_it->e_target < in3_it->e_source) {
          ++out1_it;
      } else if (in3_it->e_source < out1_it->e_target) {
          ++in3_it;
          if (in3_it == in3end) break;
          if (blockstart != out1end) {
              if (in3_it->e_source == last_e2) {
                  out1_it = blockstart;
              } else {
                  blockstart = out1end;
              }
          }
      } else {
          if (blockstart == out1end) {
              blockstart = out1_it;
          }
          last_e2 = out1_it->e_target;
          *result_it = { .rat12 = out1_it->rat_out, .e2 = last_e2, .rat23 = in3_it->rat_in };
          if (debug) cout << "      angle: " << e2label[e1] << " " << rat2label[result_it->rat12] << " "
                  << e2label[last_e2] << " " << rat2label[result_it->rat23] << " " << e2label[e3] << endl;
          ++result_it;
          ++out1_it;
      }
  }
  // shorten result to actual length and return it:
  result.resize(result_it - result.begin());
  return result;
}

#endif
