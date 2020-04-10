/** Handling of angles.
 *
 * An angle is an indirect connection between two entities e1, e3 via a middle entity e2.
 * It is specified as a triple {rat12, e2, rat23} where rat12, rat23 are two relationship or action types.
 *
 * @author Jobst Heitzig, Potsdam Institute for Climate Impact Research, heitzig@pik-potsdam.de
 * @date Mar 30, 2020
 *
 * @file
 */

// TODO: add influence of legs!

#include <assert.h>
#include <iostream>

#include "global_variables.h"
#include "probability.h"

#include "debugging.h"
#include "event.h"
#include "angle.h"

using namespace std;

void add_or_delete_angle (
        event_class ec_angle, ///< [in] event class: EC_EST adds an angle, EC_TERM deletes one
        entity e1, entity_type et1, relationship_or_action_type rat12,
        entity e2, entity_type et2, relationship_or_action_type rat23,
        entity e3, entity_type et3)
{
    if (debug) cout << "    " << ec2label[ec_angle] << " angle \"" << e2label[e1] << " " << rat2label[rat12] << " "
            << e2label[e2] << " " << rat2label[rat23] << " " << e2label[e3] << "\"" << endl;
    if ((e1 != e2) && (e2 != e3) && (e3 != e1)) {
        n_angles += (ec_angle == EC_EST) ? 1 : -1;
    }
    for (auto& rat13 : ets2relations[{ et1, et3 }]) {
        event_class ec13 = (e2outs[e1].count({ .rat_out = rat13, .e_other = e3 }) == 0) ? EC_EST : EC_TERM;
        event ev = { .ec = ec13, .e1 = e1, .rat13 = rat13, .e3 = e3 };
        event_type evt = { .ec = ec13, .et1 = et1, .rat13 = rat13, .et3 = et3 };
        if (debug) cout << "     possibly updating event: " << ec2label[ec13] <<  " \"" << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << "\"" << endl;
        if (evt2base_probunit.count(evt) > 0) {
            if (debug) cout << "      event type " << evt << " has a base success prob. of "
                    << probunit2probability(evt2base_probunit.at(evt), evt2left_tail.at(evt), evt2right_tail.at(evt)) << endl;
            influence_type inflt = {
                    .evt = evt,
                    .at = { .rat12 = rat12, .et2 = et2, .rat23 = rat23 }
            };
            auto dar = _inflt2attempt_rate[INFLT(inflt)];
            auto dsl = _inflt2delta_probunit[INFLT(inflt)];
            if (COUNT_ALL_ANGLES || (dar != 0.0) || (dsl != 0.0)) {
                if (debug) cout << "       angle may influence attempt or success" << endl;
                if (ec_angle == EC_EST) { // angle is added:
                    if (ev2data.count(ev) == 0) {
                        if (debug) cout << "        event will be scheduled newly" << endl;
                        auto evd_ = &ev2data[ev];  // generates a new event_data object
                        evd_->n_angles = 1;
                        evd_->attempt_rate = dar;
                        evd_->success_probunits = dsl;
                        schedule_event(ev, evd_, evt2left_tail.at(evt), evt2right_tail.at(evt));
                    } else {
                        if (debug) cout << "        event will be rescheduled" << endl;
                        auto evd_ = &(ev2data.at(ev));
                        evd_->n_angles += 1;
                        evd_->attempt_rate += dar;
                        evd_->success_probunits += dsl;
                        reschedule_event(ev, evd_, evt2left_tail.at(evt), evt2right_tail.at(evt));
                    }
                } else { // angle is removed:
                    auto evd_ = &(ev2data.at(ev));
                    assert (evd_->n_angles > 0);
                    evd_->n_angles -= 1;
                    evd_->attempt_rate = max(0.0, evd_->attempt_rate - dar);
                    evd_->success_probunits -= dsl;
                    if ((ec13 != EC_TERM) && (evd_->n_angles == 0)) { // only spontaneous non-termination event is left:
                        if (debug) cout << "        last angle was removed, so event will be removed" << endl;
                        // remove specific event:
                        remove_event(ev, evd_);
                    } else {
                        if (debug) cout << "        event will be rescheduled" << endl;
                        reschedule_event(ev, evd_, evt2left_tail.at(evt), evt2right_tail.at(evt));
                    }
                }
            } else {
                if (debug) cout << "       but angle may not influence attempt or success" << endl;
            }
        }
    }
}

// adapted from set_intersection template:
angles leg_intersection (entity e1, outleg_set& out1, inleg_set& in3, entity e3)
{
  auto out1_it = out1.begin();
  auto in3_it = in3.begin();
  angles result(max(out1.size(), in3.size()) * n_rats * n_rats); // allocate max. size of result
  auto result_it = result.begin();
  /** algorithm:
   *
   * the two sequences are sorted by e2.
   * put blockstart = out1.end().
   * repeat:
   *   if out1.e2 < in3.e2, advance out1.
   *   else if out1.e2 > in3.e2:
   *     advance in3.
   *     if blockstart != out1.end():
   *       if in3.e2 == previous in3.e2, rewind out2 to blockstart
   *       else put blockstart = out1.end().
   *   else out1.e2 == in3.e2:
   *     if blockstart == out1.end(), remember out1 position as blockstart
   *     store found angle
   *     advance out1
   *
   */
  auto out1end = out1.end(), blockstart = out1end;
  auto in3end = in3.end();
  entity last_e2 = out1_it->e_other;
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
              if (in3_it->e_other == last_e2) {
                  out1_it = blockstart;
              } else {
                  break;
              }
          }
      }
      if (debug) cout << "         checking: " << rat2label[out1_it->rat_out] << " " << e2label[out1_it->e_other] << ", "
              << e2label[in3_it->e_other] << " " << rat2label[in3_it->rat_in] << endl;
      if (out1_it->e_other < in3_it->e_other) {
          ++out1_it;
      } else if (in3_it->e_other < out1_it->e_other) {
          ++in3_it;
          if (in3_it == in3end) break;
          if (blockstart != out1end) {
              if (in3_it->e_other == last_e2) {
                  out1_it = blockstart;
              } else {
                  blockstart = out1end;
              }
          }
      } else {
          if (blockstart == out1end) {
              blockstart = out1_it;
          }
          last_e2 = out1_it->e_other;
          *result_it = { .rat12 = out1_it->rat_out, .e2 = last_e2, .rat23 = in3_it->rat_in };
          if (debug) cout << "      angle: " << e2label[e1] << " " << rat2label[result_it->rat12] << " "
                  << e2label[last_e2] << " " << rat2label[result_it->rat23] << " " << e2label[e3] << endl;
          ++result_it;
          ++out1_it;
      }
  }
  result.resize(result_it - result.begin()); // shorten result to actual length
  return result;
}



