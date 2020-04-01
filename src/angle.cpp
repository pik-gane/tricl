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

#include <iostream>

#include "global_variables.h"
#include "probability.h"

#include "event.h"
#include "angle.h"

using namespace std;

void add_or_delete_angle (
        event_class ec, ///< [in] event class
        entity e1, entity_type et1, relationship_or_action_type rat12,
        entity e2, entity_type et2, relationship_or_action_type rat23,
        entity e3, entity_type et3)
{
    if (debug) cout << "    " << ec2label[ec] << " angle \"" << e2label[e1] << " " << rat2label[rat12] << " "
            << e2label[e2] << " " << rat2label[rat23] << " " << e2label[e3] << "\"" << endl;
    for (auto& rat13 : ets2relations[{ et1, et3 }]) {
        event_class ec13 = (e2outs[e1].count({ .e = e3, .r = rat13 }) == 0) ? EC_EST : EC_TERM;
        event ev = { .ec = ec13, .e1 = e1, .rat13 = rat13, .e3 = e3 };

        // FIXME: deal with case where event is not scheduled yet!!!

        event_type evt = { .ec = ec13, .et1 = et1, .rat13 = rat13, .et3 = et3 };
        if (debug) cout << "     possibly updating event: " << ec2label[ec13] <<  " \"" << e2label[e1] << " " << rat2label[rat13] << " " << e2label[e3] << "\"" << endl;
        if (evt2base_probunit.count(evt) > 0) {
            if (debug) cout << "      event type " << evt << " has a base success prob. of " << probunit2probability(evt2base_probunit[evt]) << endl;
            influence_type inflt = {
                    .evt = evt,
                    .at = { .rat12 = rat12, .et2 = et2, .rat23 = rat23 }
            };
            auto dar = _inflt2attempt_rate[INFLT(inflt)], dsl = _inflt2delta_probunit[INFLT(inflt)];
            if (dar != 0.0 || dsl != 0.0) {
                if (debug) cout << "      angle may influence attempt or success" << endl;
                if (ec == EC_EST) { // angle is added:
                    if (ev2data.count(ev) == 0) {
                        if (debug) cout << "       event will be scheduled newly" << endl;
                        auto evd_ = &ev2data[ev];  // generates a new event_data object
                        evd_->n_angles = 1;
                        evd_->attempt_rate = dar;
                        evd_->success_probunit = dsl;
                        schedule_event(ev, evd_);
                    } else {
                        if (debug) cout << "       event was scheduled already, will be rescheduled" << endl;
                        auto evd_ = &ev2data.at(ev);
                        evd_->n_angles += 1;
                        evd_->attempt_rate += dar;
                        evd_->success_probunit += dsl;
//                        auto& told = evd_->t;
//                        if (told > -INFINITY) t2be.erase(told);
                        reschedule_event(ev, evd_);
                    }
                } else { // angle is removed:
                    auto evd_ = &ev2data.at(ev);
                    evd_->n_angles -= 1;
                    evd_->attempt_rate -= dar;
                    evd_->success_probunit -= dsl;
                    if ((ec13 != EC_TERM) && (evd_->n_angles == 0)) { // only spontaneous non-termination event is left:
                        // remove specific event:
                        remove_event(ev, evd_); // FIXME: sometimes event was not scheduled before...
                    } else {
//                        auto& told = evd_->t;
//                        if (told > -INFINITY) t2be.erase(told);
                        schedule_event(ev, evd_);
                    }
                }
            } else {
                if (debug) cout << "       but angle may not influence attempt or success" << endl;
            }
        }
    }
}

// adapted from set_intersection template:
angles leg_intersection (entity e1, leg_set& out, leg_set& in, entity e3)
{
  auto out_it = out.begin(), in_it = in.begin();
  angles result(out.size() + in.size()); // allocate max. size of result
  auto result_it = result.begin();
  while ((out_it != out.end()) && (in_it != in.end()))
  {
      if (out_it->e < in_it->e) {
          ++out_it;
      } else if (in_it->e < out_it->e) {
          ++in_it;
      } else {
          *result_it = { .rat12 = out_it->r, .e2 = out_it->e, .rat23 = in_it->r };
          if (debug) cout << "      angle: " << e2label[e1] << " " << rat2label[out_it->r] << " " << e2label[out_it->e] << " " << rat2label[in_it->r] << " " << e2label[e3] << endl;
          ++result_it; ++out_it; ++in_it;
      }
  }
  result.resize(result_it - result.begin()); // shorten result to actual length
  return result;
}



