/*
 * config.cpp
 *
 *  Created on: Mar 27, 2020
 *      Author: heitzig
 */

#include <iostream>

#include "data_model.h"
#include "global_variables.h"
#include "yaml-cpp/yaml.h"

using namespace std;

void read_config () {
    YAML::Node c = YAML::LoadFile("parameters.yaml"), n;

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
    // seed = n["seed"] ? n["seed"].as<???>() : 0;


//    for (auto& etspec : config["entity types"]){
//      std::cout << "entity type:" << etspec << endl;
//    }

    // const std::string username = config["username"].as<std::string>();
}





