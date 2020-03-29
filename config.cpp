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

void read_config () {
    cout << "READING CONFIG file parameters.yaml" << endl;
    YAML::Node c = YAML::LoadFile("parameters.yaml"), n, n1, n2, n3;

    // metadata (mandatory):
    n = c["metadata"];
    // meta_name = n["name"].as<string>();
    // meta_description = n["description"].as<string>();
    // meta_author = n["author"].as<string>();
    // meta_date = n["date"].as<string>();
    // meta_email = n["email"].as<string>();

    // files (mandatory):
    n = c["files"];
    gexf_filename = n["gexf"].as<string>();
    // log_filename = n["log"].as<string>();

    // limits (at least one):
    n = c["limits"];
    max_t = n["t"] ? n["t"].as<timepoint>() : INFINITY;
    // max_n_events = n["events"] ? n["events"].as<int>() : INFINITY;
    // max_wall_time = n["wall"] ? n["wall"].as<float>() : INFINITY;

    // options:
    n = c["options"];
    quiet = n["quiet"] ? n["quiet"].as<bool>() : false;
    verbose = n["verbose"] ? n["verbose"].as<bool>() : false;
    debug = n["debug"] ? n["debug"].as<bool>() : false;
    // seed = n["seed"]     ? n["seed"].as<???>() : 0;

    // entities:
    map<string, entity_type> label2et = {};
    map<string, entity> label2e = {};
    entity_type et = 1;
    entity e = 0;
    n1 = c["entities"];
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

    // action types:
    n1 = c["action types"];
    for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
        cout << "sorry, actions not supported yet!" << endl;
        assert (false);
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

    // TODO: past actions

    // named initial links:
    n = c["initial links"]["named"];
    if (n) {
        cout << "initial links" << endl;
        for (YAML::const_iterator it = n.begin(); it != n.end(); ++it) {
            string elabel1 = (*it)[0].as<string>(), ratlabel = (*it)[1].as<string>(), elabel3 = (*it)[2].as<string>();
            auto rat = label2rat[ratlabel];
            auto e1 = label2e[elabel1], e3 = label2e[elabel3];
            initial_links.insert({ e1, rat, e3 });
            cout << " " << elabel1 << " " << ratlabel << " " << elabel3 << endl;
            if (r_is_action_type[rat]) {
                float impact = (*it)[3].as<float>(); // TODO: use!!
            }
        }
    }

    // random initial links:
    n1 = c["initial links"]["random"];
}
/*
past actions:  # [negative time, entity, action type, entity]
  - [-2, Alice, posts, bullshit]

initial links:  # in "filter" mode, initial link types to be estimated need to be specified with value "estimate"
  named:
    - [Alice, is friends with, Bob]  # [entity, relationship type, entity]
    - [Celia, posts, bullshit, 10.0]  # [entity, action type, entity, impact]
  random:
    [user, is friends with, user]:  # [entity type, relationship type, entity type]
      blocks: 30  # -> block model
      within: 0.9  # link prob. in a block
      between: 0.1  # link prob. between blocks
    [user, follows, user]:
      dimension: 2  # -> random geometric graph
      decay: 3.5  # rate of exponential decay of link prob. with eucl. dist.
    # or [user, follows, user]: estimate
    [user, posts, message]:  # [entity type, relationship type, entity type] -> need a "rate"
      density: 0.01  # -> random graph
      rate: 2  # past occurrence rate of action (gives initial impacts)


dynamics:  # all values default to 0.0. in "filter" mode, parameters to be estimated need to be specified with value "estimate"
  [user, is friends with, user]:  # relations need "establish" and/or "terminate"
    establish:
      attempt:
        basic: 0.001
        [~, is friends with, user, is friends with, ~]: 1.0  # angle
      success:  # defaults to 0.0
        tails: [1.0, 2.0]
        basic: 1.0
    terminate:
      attempt:
        [~, is friends with, user]: 0.1  # outgoing leg of source entity
        [user, is friends with, ~]: 0.1  # incoming leg of destination entity
      success:
        tails: 3.0  # sets both left and right tail
  [user, posts, message]:  # actions have directly "attempt" and "success"
    attempt: ~  # as above
    success: ~  # as above


*/





