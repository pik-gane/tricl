/*
 * config.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

// for release mode, uncomment this:
//#define NDEBUG
#include <assert.h>
//

#include <iostream>

#include "data_model.h"
#include "global_variables.h"
#include "yaml-cpp/yaml.h"

using namespace std;

#include "data_model.h"
#include "global_variables.h"

bool verbose = false, quiet = false, debug = false;
timepoint max_t = 0.0;
string gexf_filename = "";
unordered_map<entity_type, label> et2label = {};
unordered_map<entity_type, entity> et2n = {};
unordered_map<relationship_or_action_type, label> rat2label = {};
unordered_map<relationship_or_action_type, bool> r_is_action_type = {};
unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv = {};
unordered_map<entity, entity_type> e2et = {};
unordered_map<entity, label> e2label = {};
set<link> initial_links = {};
unordered_map<entity_type, int> et2n_blocks = {};
unordered_map<link_type, probability> lt2initial_prob_within = {};
unordered_map<link_type, probability> lt2initial_prob_between = {};
unordered_map<entity_type, int> et2dim = {};
unordered_map<link_type, probability> lt2spatial_decay = {};
unordered_map<influence_type, rate> inflt2attempt_rate = {};
unordered_map<event_type, probunit> evt2base_probunit = {};
unordered_map<influence_type, probunit> inflt2delta_probunit = {};

void quit() {
    assert (false);
}

void read_config () {
    cout << "READING CONFIG file " << config_filename << endl;
    YAML::Node c = YAML::LoadFile(config_filename), n, n1, n2, n3;

    // metadata (mandatory):
    n = c["metadata"];
    // meta_name = n["name"].as<string>();
    // meta_description = n["description"].as<string>();
    // meta_author = n["author"].as<string>();
    // meta_date = n["date"].as<string>();
    // meta_email = n["email"].as<string>();

    // files (mandatory):
    n = c["files"];
    assert (n.IsMap());
    gexf_filename = n["gexf"].as<string>();
    // log_filename = n["log"].as<string>();

    // limits (at least one):
    n = c["limits"];
    assert (n.IsMap());
    max_t = n["t"] ? n["t"].as<timepoint>() : INFINITY;
    // max_n_events = n["events"] ? n["events"].as<int>() : INFINITY;
    // max_wall_time = n["wall"] ? n["wall"].as<float>() : INFINITY;

    // options:
    n = c["options"];
    if (n) {
        assert (n.IsMap());
        quiet = n["quiet"] ? n["quiet"].as<bool>() : false;
        verbose = n["verbose"] ? n["verbose"].as<bool>() : false;
        debug = n["debug"] ? n["debug"].as<bool>() : false;
        // seed = n["seed"]     ? n["seed"].as<???>() : 0;
    } else { // defaults:
        quiet = verbose = false;
        // seed = 0;
    }

    // entities:
    map<string, entity_type> label2et = {};
    map<string, entity> label2e = {};
    entity_type et = 1;
    entity e = 0;
    n1 = c["entities"];
    assert (n1.IsMap());
    for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
        auto etlabel = (it1->first).as<string>();
        cout << " entity type " << et << ": " << etlabel << endl;
        label2et[etlabel] = et;
        et2label[et] = etlabel;
        n2 = it1->second;
        if (n2.IsMap()) {
            et2n[et] = n2["n"].as<entity>();
            n3 = n2["labels"];
            if (n3) {
                for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                    auto elabel = it3->as<string>();
                    cout << "  entity " << e << ": " << elabel << endl;
                    label2e[elabel] = e;
                    e2label[e] = elabel;
                    e2et[e] = et;
                    e++;
                }
            }
        } else {
            et2n[et] = n2.as<entity>();
        }
        et++;
    }

    map<string, relationship_or_action_type> label2rat = { {"=", RT_ID} };
    rat2label = { {RT_ID, "="} };
    r_is_action_type = { {RT_ID, false} };
    rat2inv = { {RT_ID, RT_ID} };
    relationship_or_action_type nextrat = RT_ID + 1, rat1, rat2;

    // relationship types:
    n1 = c["relationship types"];
    if (n1) {
        assert (n1.IsMap());
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto ratlabel = (it1->first).as<string>();
            label2rat[ratlabel] = nextrat;
            rat2label[nextrat] = ratlabel;
            r_is_action_type[nextrat] = false;
            nextrat++;
        }
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto ratlabel = (it1->first).as<string>();
            assert (ratlabel != "symmetric");
            rat1 = label2rat[ratlabel];
            n2 = it1->second;
            if (n2.IsNull()) {
                cout << " relationship type " << rat1 << ": " << ratlabel << endl;
            } else {
                auto ratlabel2 = n2.as<string>();
                if (ratlabel2 == "symmetric") ratlabel2 = ratlabel;
                if (label2rat.count(ratlabel2)) {
                    rat2 = label2rat[ratlabel2];
                    cout << " relationship type " << rat1 << ": " << ratlabel << " (inverse: " << ratlabel2 << ")" << endl;
                } else {
                    rat2 = nextrat;
                    label2rat[ratlabel2] = rat2;
                    rat2label[rat2] = ratlabel2;
                    r_is_action_type[nextrat] = false;
                    nextrat++;
                    cout << " relationship type " << rat1 << ": " << ratlabel << " (inverse: " << ratlabel2 << ")" << endl;
                    cout << " relationship type " << rat2 << ": " << ratlabel2 << " (inverse: " << ratlabel << ")" << endl;
                }
                rat2inv[rat1] = rat2;
                rat2inv[rat2] = rat1;
            }
        }
    }

    // action types:
    n1 = c["action types"];
    if (n1) {
        assert (n1.IsMap());
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            cout << "sorry, actions not supported yet!" << endl;
            quit();
            auto ratlabel = (it1->first).as<string>();
            label2rat[ratlabel] = nextrat;
            rat2label[nextrat] = ratlabel;
            r_is_action_type[nextrat] = true;
            nextrat++;
        }
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto ratlabel = (it1->first).as<string>();
            assert (ratlabel != "symmetric");
            rat1 = label2rat[ratlabel];
            n2 = it1->second;
            if (n2.IsNull()) {
                cout << " action type " << rat1 << ": " << ratlabel << endl;
            } else {
                auto ratlabel2 = n2.as<string>();
                if (ratlabel2 == "symmetric") ratlabel2 = ratlabel;
                if (label2rat.count(ratlabel2)) {
                    rat2 = label2rat[ratlabel2];
                    cout << " action type " << rat1 << ": " << ratlabel << " (inverse: " << ratlabel2 << ")" << endl;
                } else {
                    rat2 = nextrat;
                    label2rat[ratlabel2] = rat2;
                    rat2label[rat2] = ratlabel2;
                    r_is_action_type[nextrat] = true;
                    nextrat++;
                    cout << " action type " << rat1 << ": " << ratlabel << " (inverse: " << ratlabel2 << ")" << endl;
                    cout << " action type " << rat2 << ": " << ratlabel2 << " (inverse: " << ratlabel << ")" << endl;
                }
                rat2inv[rat1] = rat2;
                rat2inv[rat2] = rat1;
            }
        }
    }

    // TODO: verify ets, es and rats mentioned later exist!

    // TODO: past actions
    /*
    past actions:  # [negative time, entity, action type, entity]
      - [-2, Alice, posts, bullshit]
    */

    // named initial links:
    n = c["initial links"]["named"];
    if (n) {
        cout << "initial links" << endl;
        for (YAML::const_iterator it = n.begin(); it != n.end(); ++it) {
            assert (it->IsSequence());
            string elabel1 = (*it)[0].as<string>(),
                    ratlabel = (*it)[1].as<string>(),
                    elabel3 = (*it)[2].as<string>();
            auto rat = label2rat.at(ratlabel);
            auto e1 = label2e.at(elabel1), e3 = label2e.at(elabel3);
            initial_links.insert({ e1, rat, e3 });
            cout << " " << elabel1 << " " << ratlabel << " " << elabel3 << endl;
            if (r_is_action_type[rat]) {
//                float impact = (*it)[3].as<float>(); // TODO: use!!
            }
        }
    }

    // random initial links:
    n1 = c["initial links"]["random"];
    if (n1) {
        assert (n1.IsMap());
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto lt = it1->first;
            assert (lt.IsSequence());
            auto et1 = label2et.at(lt[0].as<string>()),
                    et3 = label2et.at(lt[2].as<string>());
            auto rat = label2rat.at(lt[1].as<string>());
            auto spec = it1->second;
            assert (spec.IsMap());
            if (spec["density"]) { // Erdös-Renyi random graph, treated as block model with one block
                auto pw = spec["density"].as<probability>();
                assert ((0.0 <= pw) && (pw <= 1.0));
//                auto rate = spec["rate"] ? spec["rate"].as<float>() : 1.0; // TODO: for action types
                et2n_blocks[et1] = et2n_blocks[et3] = 1;
                lt2initial_prob_within[{et1, rat, et3}] = pw;
                lt2initial_prob_between[{et1, rat, et3}] = 0.0;
            } else if (spec["blocks"]) {
                if (et1 == et3) { // symmetric block model
                    auto n = spec["blocks"].as<int>();
                    // TODO: allow list of "sizes"
                    auto pw = spec["within"] ? spec["within"].as<probability>() : 1.0;
                    auto pb = spec["between"] ? spec["between"].as<probability>() : 0.0;
                    assert ((n > 0) && (0.0 <= pw) && (pw <= 1.0) && (0.0 <= pb) && (pb <= 1.0));
                    et2n_blocks[et1] = n;
                    lt2initial_prob_within[{et1, rat, et3}] = pw;
                    lt2initial_prob_between[{et1, rat, et3}] = pb;
                } else { // asymmetric block model
                    // TODO
                    cout << "sorry, block model for asymmetric relationship or action types not supported yet!" << endl;
                    quit();
                }
            } else if (spec["dimension"]) {
                // TODO: make sure no conflicts between several spatial models for same entity types!
                auto dim = spec["dimension"].as<int>();
                // TODO: allow list of "widths"
                auto dec = spec["decay"] ? spec["decay"].as<float>() : 1.0;
                assert ((dim > 0) && (dec > 0.0));
                et2dim[et1] = et2dim[et3] = dim;
                lt2spatial_decay[{et1, rat, et3}] = dec;
            }
        }
    }
    // random initial links:
    n1 = c["dynamics"];
    assert (n1.IsMap());
    for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
        auto lt = it1->first;
        assert (lt.IsSequence());
        auto et1 = label2et.at(lt[0].as<string>()),
                et3 = label2et.at(lt[2].as<string>());
        auto rat13 = label2rat.at(lt[1].as<string>());
        auto spec = it1->second;
        assert (spec.IsMap());
        if (r_is_action_type[rat13]) {
            // TODO: attempt, success
        } else {
            auto n2 = spec["establish"];
            auto ec = EC_EST;
            evt2base_probunit[{ ec, et1, rat13, et3 }] = 0.0;
            if (n2) {
                assert (n2.IsMap());
                auto n3 = n2["attempt"];
                if (n3) {
                    assert (n3.IsMap());
                    for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                        auto cause = it3->first;
                        auto ar = it3->second.as<rate>();
                        assert (ar >= 0.0);
                        entity_type et2;
                        relationship_or_action_type rat12, rat23;
                        if (cause.IsSequence()) {
                            if (cause.size() == 5) { // angle
                                assert (cause[0].IsNull() && cause[4].IsNull());
                                rat12 = label2rat.at(cause[1].as<string>());
                                et2 = label2et.at(cause[2].as<string>());
                                rat23 = label2rat.at(cause[3].as<string>());
                            } else {
                                // TODO later
                                cout << "sorry, legs cannot attempt establishment yet" << endl;
                                quit();
                                assert (cause.size() == 3);
                                if (cause[0].IsNull()) { // outgoing leg
                                    rat12 = label2rat.at(cause[1].as<string>());
                                    et2 = label2et.at(cause[2].as<string>());
                                    rat23 = NO_RAT;
                                } else { // incoming leg
                                    rat12 = NO_RAT;
                                    et2 = label2et.at(cause[0].as<string>());
                                    rat23 = label2rat.at(cause[1].as<string>());
                                    assert (cause[2].IsNull());
                                }
                            }
                            inflt2attempt_rate[{ { ec, et1, rat13, et3 }, { rat12, et2, rat23 } }] = ar;
                            inflt2delta_probunit[{ { ec, et1, rat13, et3 }, { rat12, et2, rat23 } }] = 0.0;
                        } else { // basic
                            assert (cause.as<string>() == "basic");
                            inflt2attempt_rate[{ { ec, et1, rat13, et3 }, NO_ANGLE }] = ar;
                        }
                    }
                }
                n3 = n2["success"];
                if (n3) {
                    // TODO: tails!
                    assert (n3.IsMap());
                    for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                        auto cause = it3->first;
                        auto pu = it3->second.as<probunit>();
                        entity_type et2;
                        relationship_or_action_type rat12, rat23;
                        if (cause.IsSequence()) {
                            if (cause.size() == 5) { // angle
                                assert (cause[0].IsNull() && cause[4].IsNull());
                                rat12 = label2rat.at(cause[1].as<string>());
                                et2 = label2et.at(cause[2].as<string>());
                                rat23 = label2rat.at(cause[3].as<string>());
                            } else {
                                // TODO later
                                cout << "sorry, legs cannot influence establishment yet" << endl;
                                quit();
                                assert (cause.size() == 3);
                                if (cause[0].IsNull()) { // outgoing leg
                                    rat12 = label2rat.at(cause[1].as<string>());
                                    et2 = label2et.at(cause[2].as<string>());
                                    rat23 = NO_RAT;
                                } else { // incoming leg
                                    rat12 = NO_RAT;
                                    et2 = label2et.at(cause[0].as<string>());
                                    rat23 = label2rat.at(cause[1].as<string>());
                                    assert (cause[2].IsNull());
                                }
                            }
                            inflt2delta_probunit[{ { ec, et1, rat13, et3 }, { rat12, et2, rat23 } }] = pu;
                        } else { // basic
                            assert (cause.as<string>() == "basic");
                            evt2base_probunit[{ ec, et1, rat13, et3 }] = pu;
                        }
                    }
                }
            }
            n2 = spec["terminate"];
            ec = EC_TERM;
            evt2base_probunit[{ ec, et1, rat13, et3 }] = 0.0;
            if (n2) {
                assert (n2.IsMap());
                auto n3 = n2["attempt"];
                if (n3) {
                    assert (n3.IsMap());
                    for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                        auto cause = it3->first;
                        auto ar = it3->second.as<rate>();
                        assert (ar >= 0.0);
                        entity_type et2;
                        relationship_or_action_type rat12, rat23;
                        if (cause.IsSequence()) {
                            if (cause.size() == 5) { // angle
                                assert (cause[0].IsNull() && cause[4].IsNull());
                                rat12 = label2rat.at(cause[1].as<string>());
                                et2 = label2et.at(cause[2].as<string>());
                                rat23 = label2rat.at(cause[3].as<string>());
                            } else {
                                assert (cause.size() == 3);
                                if (cause[0].IsNull()) { // outgoing leg
                                    rat12 = label2rat.at(cause[1].as<string>());
                                    et2 = label2et.at(cause[2].as<string>());
                                    rat23 = NO_RAT;
                                } else { // incoming leg
                                    rat12 = NO_RAT;
                                    et2 = label2et.at(cause[0].as<string>());
                                    rat23 = label2rat.at(cause[1].as<string>());
                                    assert (cause[2].IsNull());
                                }
                            }
                            inflt2attempt_rate[{ { ec, et1, rat13, et3 }, { rat12, et2, rat23 } }] = ar;
                            inflt2delta_probunit[{ { ec, et1, rat13, et3 }, { rat12, et2, rat23 } }] = 0.0;
                        } else { // basic
                            assert (cause.as<string>() == "basic");
                            inflt2attempt_rate[{ { ec, et1, rat13, et3 }, NO_ANGLE }] = ar;
                        }
                    }
                }
                n3 = n2["success"];
                if (n3) {
                    // TODO: tails!
                    assert (n3.IsMap());
                    for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                        auto cause = it3->first;
                        auto pu = it3->second.as<probunit>();
                        entity_type et2;
                        relationship_or_action_type rat12, rat23;
                        if (cause.IsSequence()) {
                            if (cause.size() == 5) { // angle
                                assert (cause[0].IsNull() && cause[4].IsNull());
                                rat12 = label2rat.at(cause[1].as<string>());
                                et2 = label2et.at(cause[2].as<string>());
                                rat23 = label2rat.at(cause[3].as<string>());
                            } else {
                                assert (cause.size() == 3);
                                if (cause[0].IsNull()) { // outgoing leg
                                    rat12 = label2rat.at(cause[1].as<string>());
                                    et2 = label2et.at(cause[2].as<string>());
                                    rat23 = NO_RAT;
                                } else { // incoming leg
                                    rat12 = NO_RAT;
                                    et2 = label2et.at(cause[0].as<string>());
                                    rat23 = label2rat.at(cause[1].as<string>());
                                    assert (cause[2].IsNull());
                                }
                            }
                            inflt2delta_probunit[{ { ec, et1, rat13, et3 }, { rat12, et2, rat23 } }] = pu;
                        } else { // basic
                            assert (cause.as<string>() == "basic");
                            evt2base_probunit[{ ec, et1, rat13, et3 }] = pu;
                        }
                    }
                }
            }
        }
    }
}

