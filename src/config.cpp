/** Handling of configuration files.
 *
 *  \file
 *
 *  Configuration files are in YAML.
 *  See the repository's \ref README.md for a documentation of all options.
 *
 */

#include <limits.h>
#include <math.h>
#include <iostream>
#include "yaml-cpp/yaml.h"
#include "3rdparty/tinyexpr.h"   // for handling of expressions
#include "3rdparty/cxxopts.hpp"  // for handling of command line options

#include "global_variables.h"
#include "entity.h"
#include "io.h"

cxxopts::Options options("tricl", "a generic network-based social simulation model");  ///< Holds all command line options

// scalar parameters and their default values:
unordered_map<relationship_or_action_type, string> gexf_filename = {};
string diagram_fileprefix = "", gexf_default_filename = "";
bool verbose = false, quiet = false, debug = false;
timepoint max_t = 0.0;
long int max_n_events = LONG_MAX;
unsigned seed = 0;

// maps and sets of parameters with some defaults:
unordered_map<entity_type, label> et2label = {};
unordered_map<string, entity_type> label2et = {};
unordered_map<entity_type, entity> et2n = {};
unordered_map<relationship_or_action_type, label> rat2label = { {RT_ID, "="} };
unordered_map<string, relationship_or_action_type> label2rat = { {"=", RT_ID} };
unordered_map<relationship_or_action_type, bool> r_is_action_type = {};
unordered_map<relationship_or_action_type, relationship_or_action_type> rat2inv = {};
unordered_map<entity, label> e2label = {};
unordered_map<string, entity> label2e = {};
set<link> initial_links = {};
unordered_map<entity_type, int> et2n_blocks = {};
unordered_map<link_type, probability> lt2initial_prob_within = {};
unordered_map<link_type, probability> lt2initial_prob_between = {};
unordered_map<entity_type, int> et2dim = {};
unordered_map<link_type, probability> lt2spatial_decay = {};
unordered_map<event_type, rate> evt2base_attempt_rate = {};
unordered_map<influence_type, rate> inflt2attempt_rate = {};
unordered_map<event_type, double> evt2left_tail = {}, evt2right_tail = {};
unordered_map<event_type, probunits> evt2base_probunits = {};
unordered_map<influence_type, probunits> inflt2delta_probunits = {};
unordered_map<entity_type, double> et2gexf_size = {}, et2gexf_a = {};
unordered_map<entity_type, string> et2gexf_shape = {};
unordered_map<entity_type, int> et2gexf_r = {}, et2gexf_g = {}, et2gexf_b = {};
unordered_map<relationship_or_action_type, double> rat2gexf_thickness = {}, rat2gexf_a = {};
unordered_map<relationship_or_action_type, string> rat2gexf_shape = {};
unordered_map<relationship_or_action_type, int> rat2gexf_r = {}, rat2gexf_g = {}, rat2gexf_b = {};

// further default values:
#define TAIL_DEFAULT 1.0  ///< Default value for tail indices. TODO: change to 0.0 to increase default performance?

// internal data:
entity max_e = 0;

// stuff for parsing expressions:
#define MAX_N_TE_VARS 1000           // maximum number of metaparameters
string te_syms[MAX_N_TE_VARS];       // symbols representing metaparameters in expressions
double te_vals[MAX_N_TE_VARS];       // corresponding substitution values
te_variable te_vars[MAX_N_TE_VARS];  // corresponding variable objects
int n_te_vars = 0;                   // total no. of these

// convert a string expression into a double value:
double parse_double (string expr)
{
    return te_eval(te_compile(expr.c_str(), te_vars, n_te_vars, 0));
}

// convert a string expression into an integer value:
double parse_int (string expr)
{
    double dblv = parse_double(expr);
    return (dblv == INFINITY) ? INT_MAX : (dblv == -INFINITY) ? INT_MIN : round(dblv);
}

// register a metaparameter and its value for use in expressions:
double register_te_var (string symbol, string expr)
{
    if (n_te_vars >= MAX_N_TE_VARS) throw "too many metaparameters";
    te_syms[n_te_vars] = symbol;
    // compute and append value to te_vals:
    double value = parse_double(expr);
    te_vals[n_te_vars] = {value};
    // register symbol and value as new subexpression to be substituted:
    te_vars[n_te_vars] = {te_syms[n_te_vars].c_str(), (const void*)(&te_vals[n_te_vars]), TE_VARIABLE};
    n_te_vars++;
    return value;
}

// read entity labels from a YAML list
void read_entity_labels (YAML::Node n, entity_type et)
{
    for (YAML::const_iterator it3 = n.begin(); it3 != n.end(); ++it3) {
        auto elabel = it3->as<string>();  // YAML values must always be converted to some type by means of as<type>()
        entity e = add_entity(et, elabel);
        cout << "  entity " << e << ": " << elabel << endl;
    }
}

/** Parse the command line options and YAML config file.
 */
void read_config (
        int argc,     ///< [in] forwarded from \ref main()
        char *argv[]  ///< [in] forwarded from \ref main()
        )
{
    // parse command line:

    if (argc<2) {
        cout << "Usage:  tricl config_file.yaml  OR  tricl config_file.yaml [options]  OR  tricl config_file.yaml --help";
        exit(1);
    }
    config_yaml_filename = argv[1];
    YAML::Node c = YAML::LoadFile(config_yaml_filename), n, n1, n2, n3;
    options.add_options()("help", "Print help for " + config_yaml_filename);

    // register command line options for all metaparameters in config file:
    n = c["metaparameters"];
    if (n) {
        if (!n.IsMap()) {
            cout << options.help() << endl;
            throw "field 'metaparameters' in config file must be a map";
        }
        for (YAML::const_iterator it = n.begin(); it != n.end(); ++it) {  // iterate through the map
            string symbol = it->first.as<string>();  // the "first" member of a YAML map entry is the item's key
            string deflt = it->second.as<string>();  // the "second" member of a YAML map entry is the item's value
            // TODO: validate symbol is both a valid command line option name and a valid tinyexpr variable name!
            // the expression specified in the config file acts as default
            // and can be overridden by the same named command-line option:
            options.add_options()(symbol, "value or expression for metaparameter " + symbol,
                    cxxopts::value<std::string>()->default_value(deflt));
        }
    }
    // read command line:
    auto cmdlineopts = options.parse(argc, argv);
    if ((cmdlineopts.count("help") >= 1)) {
        cout << options.help() << endl;
        exit(0);
    }

    // read config file:

    if (!quiet) cout << "READING CONFIG file " << config_yaml_filename << " ..." << endl;

    // initial symbols allowed in expressions:
    te_vals[0] = { INFINITY };
    te_vars[0] = {"inf", (const void*)(&te_vals[0]), TE_VARIABLE};
    te_vals[1] = { INFINITY };
    te_vars[1] = {"infty", (const void*)(&te_vals[0]), TE_VARIABLE};
    te_vals[2] = { INFINITY };
    te_vars[2] = {"infinity", (const void*)(&te_vals[0]), TE_VARIABLE};
    te_vals[3] = { std::numeric_limits<double>::epsilon() };
    te_vars[3] = {"eps", (const void*)(&te_vals[0]), TE_VARIABLE};
    te_vals[4] = { std::numeric_limits<double>::epsilon() };
    te_vars[4] = {"epsilon", (const void*)(&te_vals[0]), TE_VARIABLE};
    n_te_vars = 5;

    // NOTE: for some reason, "metaparameters" must be parsed before "options" to avoid a strange error.

    // metaparameters:
    n = c["metaparameters"];
    if (n) {
        if (!quiet) cout << " metaparameters:" << endl;
        for (YAML::const_iterator it = n.begin(); it != n.end(); ++it) {
            auto symbol = it->first.as<string>();
            auto expression = cmdlineopts[symbol].as<string>();  // this is either the expression from the command line, or if missing then the one from the config file
            double value = register_te_var(symbol, expression);
            if (!quiet) cout << "  " << symbol << " = " << expression << " = " << value << endl;
        }
    }

    // options:
    n = c["options"];
    if (n) {
        if (!n.IsMap()) throw "yaml field 'options' must be a map";
        if (n["quiet"]) quiet = n["quiet"].as<bool>();
        if (n["verbose"]) verbose = n["verbose"].as<bool>();
        if (n["debug"]) debug = n["debug"].as<bool>();
        if (n["seed"]) seed = n["seed"].as<unsigned>();
    }

    // TODO: make use of the metadata, e.g. in output
    // metadata (mandatory):
    n = c["metadata"];
    // meta_name = n["name"].as<string>();
    // meta_description = n["description"].as<string>();
    // meta_author = n["author"].as<string>();
    // meta_date = n["date"].as<string>();
    // meta_email = n["email"].as<string>();

    // files (mandatory):
    n = c["files"];
    if (!n.IsMap()) throw "yaml field 'files' must be a map";
    if (n["gexf"]) gexf_default_filename = n["gexf"].as<string>();
    // log_filename = n["log"].as<string>();
    if (n["diagram prefix"]) diagram_fileprefix = n["diagram prefix"].as<string>();

    // limits (at least one):
    n = c["limits"];
    if (!n.IsMap()) throw "yaml field 'limits' must be a map";
    if (n["t"]) {
        max_t = parse_double(n["t"].as<string>());
        if (max_t == INFINITY) throw "limit: t must be finite";
    }
    if (n["events"]) max_n_events = floor(parse_double(n["events"].as<string>()));
    // max_wall_time = n["wall"] ? n["wall"].as<double>() : INFINITY;  // TODO: support this option
    if ((max_t==INFINITY) && (max_n_events==LONG_MAX)) throw
            "must specify at least one of limits:t, limits:events";

    // entities:
    entity_type et = 1;
    n1 = c["entities"];
    if (!n1.IsMap()) throw "yaml field 'entities' must be a map";
    for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
        auto etlabel = (it1->first).as<string>();
        cout << " entity type " << et << ": " << etlabel << endl;
        label2et[etlabel] = et;
        et2label[et] = etlabel;
        n2 = it1->second;
        if (n2.IsMap()) {
            et2n[et] = parse_int(n2["n"].as<string>());
            n3 = n2["labels"] ? n2["labels"] : n2["names"] ? n2["names"] : n2["named"];
            if (n3) read_entity_labels(n3, et);
        } else if (n2.IsSequence()) {  // only a list of labels is specified
            // read and count them:
            entity last_e = max_e;
            read_entity_labels(n2, et);
            et2n[et] = max_e - last_e;
        } else {
            et2n[et] = parse_int(n2.as<string>());
        }
        // defaults:
        et++;
    }

    rat2label = { {RT_ID, "="} };
    r_is_action_type = { {RT_ID, false} };
    rat2inv = { {RT_ID, RT_ID} };
    relationship_or_action_type nextrat = RT_ID + 1, rat1, rat2;

    // relationship types:
    n1 = c["relationship types"];
    if (n1) {
        if (!n1.IsMap()) throw "yaml field 'relationship types' must be a map";
        // go through list once to collect labels:
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto ratlabel = (it1->first).as<string>();
            label2rat[ratlabel] = nextrat;
            rat2label[nextrat] = ratlabel;
            r_is_action_type[nextrat] = false;
            nextrat++;
        }
        // go through it again to collect further information:
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto ratlabel = (it1->first).as<string>();
            if (ratlabel == "symmetric") throw "'symmetric' is not a valid label";
            rat1 = label2rat[ratlabel];
            n2 = it1->second;
            string ratlabel2 = "";
            if (n2.IsNull()) {
                gexf_filename[rat1] = gexf_default_filename;
            } else if (n2.IsScalar()) {
                ratlabel2 = n2.as<string>();
                gexf_filename[rat1] = gexf_default_filename;
            } else if (n2.IsMap()) {
                ratlabel2 = n2["inverse"] ? n2["inverse"].as<string>() : "";
                gexf_filename[rat1] = (!n2["gexf"]) ? gexf_default_filename : n2["gexf"].IsNull() ? "" : n2["gexf"].as<string>();
            } else throw "field values under 'relationship types' must be strings or maps";
            if (ratlabel2 == "") {
                cout << " relationship type " << rat1 << ": " << ratlabel << endl;
            } else {
                relationship_or_action_type rat2;
                if (ratlabel2 == "symmetric") ratlabel2 = ratlabel;
                if (label2rat.count(ratlabel2)) {
                    rat2 = label2rat[ratlabel2];
                    cout << " relationship type " << rat1 << ": " << ratlabel << " (inverse: " << ratlabel2 << ")" << endl;
                } else {
                    rat2 = nextrat;
                    label2rat[ratlabel2] = rat2;
                    rat2label[rat2] = ratlabel2;
                    gexf_filename[rat2] = ""; // don't output rats that are only inverses
                    r_is_action_type[nextrat] = false;
                    nextrat++;
                    cout << " relationship type " << rat2 << ": " << ratlabel2 << " (inverse: " << ratlabel << ")" << endl;
                    cout << " relationship type " << rat1 << ": " << ratlabel << " (inverse: " << ratlabel2 << ")" << endl;
                }
                rat2inv[rat1] = rat2;
                rat2inv[rat2] = rat1;
            }
            if (gexf_filename[rat1] == "") {
                cout << "  (omitted in gexf output)" << endl;
            } else if (gexf_filename[rat1] != gexf_default_filename) {
                cout << "  written to separate file " << gexf_filename[rat1] << endl;
            }
        }
    }

    // action types:
    n1 = c["action types"];
    if (n1) {
        if (!n1.IsMap()) throw "yaml field 'action types' must be a map";
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            throw "sorry, actions not supported yet!";
            auto ratlabel = (it1->first).as<string>();
            label2rat[ratlabel] = nextrat;
            rat2label[nextrat] = ratlabel;
            r_is_action_type[nextrat] = true;
            nextrat++;
        }
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto ratlabel = (it1->first).as<string>();
            if (ratlabel == "symmetric") throw "'symmetric' is not a valid label";
            rat1 = label2rat[ratlabel];
            n2 = it1->second;
            if (n2.IsNull()) {
                cout << " action type " << rat1 << ": " << ratlabel << endl;
            } else {
                auto ratlabel2 = n2.as<string>();
                if (ratlabel2 == "symmetric") ratlabel2 = ratlabel;
                if (label2rat.count(ratlabel2) > 0) {
                    rat2 = label2rat[ratlabel2];
                    cout << " action type " << rat1 << ": " << ratlabel << " (inverse: " << ratlabel2 << ")" << endl;
                    if ((rat2inv.count(rat1) > 0) && (rat2inv.at(rat1) != rat2)) throw "conflicting specification of inverses";
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
        cout << "named initial links:" << endl;
        for (YAML::const_iterator it = n.begin(); it != n.end(); ++it) {
            if (!it->IsSequence()) throw "yaml subfield 'named' of 'initial links' must be a sequence";
            string e1label = (*it)[0].as<string>(),
                    rat13label = (*it)[1].as<string>(),
                    e3label = (*it)[2].as<string>();
            cout << " " << e1label << " " << rat13label << " " << e3label << endl;
            try {
                auto e1 = label2e.at(e1label), e3 = label2e.at(e3label);
                auto rat13 = label2rat.at(rat13label);
                initial_links.insert({ e1, rat13, e3 });
                if (r_is_action_type[rat13]) {
//                float impact = (*it)[3].as<double>(); // TODO: use!!
                }
            } catch (const std::exception&) {
                throw "some entity or the relationship or action type was not declared";
            }
        }
    }

    // random initial links:
    n1 = c["initial links"]["random"];
    if (n1) {
        cout << "random initial links:" << endl;
        if (!n1.IsMap()) throw "yaml subfield 'random' of 'initial links' must be a map";
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto lt = it1->first;
            if (!lt.IsSequence()) throw "keys in yaml map 'named' of 'initial links' must be of the form [entity type, relationship or action type, entity type]";
            cout << " " << lt[0].as<string>() << " " << lt[1].as<string>() << " " << lt[2].as<string>() << endl;
            try {
                auto et1 = label2et.at(lt[0].as<string>()),
                        et3 = label2et.at(lt[2].as<string>());
                auto rat = label2rat.at(lt[1].as<string>());
                auto spec = it1->second;
                if (!spec.IsMap()) throw "values in yaml map 'named' of 'initial links' must be maps";
                if (spec["density"] || spec["probability"]) { // Erdös-Renyi random graph, treated as block model with one block
                    auto pw = parse_double((spec["density"] ? spec["density"] : spec["probability"]).as<string>());
                    if (!((0.0 <= pw) && (pw <= 1.0))) throw "'density'/'probability' must be between 0.0 and 1.0";
    //                auto rate = spec["rate"] ? spec["rate"].as<double>() : 1.0; // TODO: for action types
                    et2n_blocks[et1] = et2n_blocks[et3] = 1;
                    lt2initial_prob_within[{et1, rat, et3}] = pw;
                    lt2initial_prob_between[{et1, rat, et3}] = 0.0;
                } else if (spec["blocks"]) {
                    if (et1 == et3) { // symmetric block model
                        int n = parse_int(spec["blocks"].as<string>());
                        // TODO: allow list of "sizes"
                        auto pw = spec["within"] ? parse_double(spec["within"].as<string>()) : 1.0;
                        auto pb = spec["between"] ? parse_double(spec["between"].as<string>()) : 0.0;
                        if (!(n > 0)) throw "'n' must be positive";
                        if (!((0.0 <= pw) && (pw <= 1.0))) throw "'within' must be between 0.0 and 1.0";
                        if (!((0.0 <= pb) && (pb <= 1.0))) throw "'between' must be between 0.0 and 1.0";
                        et2n_blocks[et1] = n;
                        lt2initial_prob_within[{et1, rat, et3}] = pw;
                        lt2initial_prob_between[{et1, rat, et3}] = pb;
                    } else { // asymmetric block model
                        // TODO
                        throw "sorry, block model for asymmetric relationship or action types not supported yet!";
                    }
                } else if (spec["dimension"]) {
                    // TODO: make sure no conflicts between several spatial models for same entity types!
                    int dim = parse_int(spec["dimension"].as<string>());
                    // TODO: allow list of "widths"
                    auto dec = spec["decay"] ? parse_double(spec["decay"].as<string>()) : 1.0;
                    if (!(dim > 0)) throw "'dimension' must be positive";
                    if (!(dec > 0.0)) throw "'decay' must be positive";
                    et2dim[et1] = et2dim[et3] = dim;
                    lt2spatial_decay[{et1, rat, et3}] = dec;
                }
            } catch (const std::exception&) {
                throw "some entity or the relationship or action type was not declared";
            }
        }
    }

    // initial links read from files:
    n1 = c["initial links"];
    if (n1) {
        if (!n1.IsMap()) throw "yaml field 'initial links' must be a map";
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto key = (it1->first).as<string>();
            auto n2 = it1->second;
            if ((key != "named") && (key != "random")) {
                auto filename = key;
                if (!quiet) cout << "reading initial links from file " << filename << " ..." << endl;
                auto ext = filename.substr(filename.find_last_of(".") + 1);
                if (ext == "csv") {
                    string et1label, rat13label, et3label;
                    entity_type et1_default = -1, et3_default = -1;
                    auto rat13 = NO_RAT;
                    if (n2["type"]) {
                        et1label = n2["type"][0].as<string>();
                        rat13label = n2["type"][1].as<string>();
                        et3label = n2["type"][2].as<string>();
                        if ((label2et.count(et1label) == 0) || (label2et.count(et3label) == 0)) throw
                                "unregistered entity type";
                        if (label2rat.count(rat13label) == 0) throw "unregistered relationship or action type";
                        et1_default = label2et.at(et1label);
                        rat13 = label2rat.at(rat13label);
                        et3_default = label2et.at(et3label);
                    } else {
                        if (n2["entity types"]) {
                            et1label = n2["entity types"][0].as<string>();
                            et3label = n2["entity types"][1].as<string>();
                            if ((label2et.count(et1label) == 0) || (label2et.count(et3label) == 0)) throw
                                    "unregistered entity type";
                            et1_default = label2et.at(et1label);
                            et3_default = label2et.at(et3label);
                        }
                        // TODO: action type
                        if (n2["relationship type"]) {
                            rat13label = n2["relationship type"].as<string>();
                            if (label2rat.count(rat13label) == 0) throw "unregistered relationship type";
                            rat13 = label2rat.at(rat13label);
                        }
                    }
                    int skip = n2["skip"] ? (parse_int(n2["skip"].as<string>())) : 0;
                    int max = n2["max"] ? (parse_int(n2["max"].as<string>())) : (INT_MAX - 1);
                    char delimiter = n2["delimiter"] ? n2["delimiter"].as<char>() : ',';
                    string prefix = n2["prefix"] ? n2["prefix"].as<string>() : ""; // TODO: distinguish source, target prefixes
                    auto cols = n2["cols"];
                    read_links_csv (filename, skip, max, delimiter,
                            cols[0].as<int>(),
                            (rat13 == NO_RAT) ? parse_int(cols[1].as<string>()) : -1,
                            (rat13 == NO_RAT) ? parse_int(cols[2].as<string>()) : parse_int(cols[1].as<string>()),
                            et1_default, rat13, et3_default,
                            prefix);
                } else {
                    throw "unsupported file extension";
                }
                if (!quiet) cout << " ...done." << endl;
            }
        }
    }

    // dynamics:
    cout << "dynamics:" << endl;
    n1 = c["dynamics"];
    if (!n1.IsMap()) throw "yaml field 'dynamics' must be a map";
    for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
        auto lt = it1->first;
        if (!lt.IsSequence()) throw "keys in yaml map 'dynamics' must be of the form [entity type, relationship or action type, entity type]";
        cout << " " << lt[0].as<string>() << " " << lt[1].as<string>() << " " << lt[2].as<string>() << endl;
        try {
            auto et1l = lt[0].as<string>(), et3l = lt[2].as<string>();
            auto et1 = label2et.at(et1l), et3 = label2et.at(et3l);
            auto rat13 = label2rat.at(lt[1].as<string>());
            auto spec = it1->second;
            if (!spec.IsMap()) throw "values in yaml map 'dynamics' must be maps";
            if (r_is_action_type[rat13]) {
                // TODO: attempt, success
            } else {
                auto n2 = spec["establish"];
                auto ec = EC_EST;
                if (n2) {
                    event_type evt = { ec, et1, rat13, et3 };
                    if (!n2.IsMap()) throw "yaml field 'establish' within 'dynamics' must be a map";
                    auto n3 = n2["attempt"];
                    if (n3) {
                        evt2left_tail[evt] = evt2right_tail[evt] = TAIL_DEFAULT;
                        evt2base_probunits[evt] = 0.0;
                        if (!n3.IsMap()) {
                            if (!n3.IsScalar()) throw "yaml field 'attempt' within 'dynamics' must be a scalar (base attempt rate) or a map";
                            evt2base_attempt_rate[evt] = parse_double(n3.as<string>());
                        } else {
                            for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                                auto cause = it3->first;
                                auto ar = parse_double(it3->second.as<string>());
                                if (!(ar >= 0.0)) throw "values in map 'attempt' must be non-negative";
                                entity_type et2;
                                relationship_or_action_type rat12, rat23;
                                if (cause.IsSequence()) {
                                    if (cause.size() == 5) { // angle
                                        if (!((cause[0].IsNull() || (cause[0].as<string>() == et1l)) && (cause[4].IsNull() || (cause[4].as<string>() == et3l)))) throw
                                                "keys in map 'attempt' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        rat12 = label2rat.at(cause[1].as<string>());
                                        et2 = label2et.at(cause[2].as<string>());
                                        rat23 = label2rat.at(cause[3].as<string>());
                                    } else {
                                        // TODO later
                                        throw "sorry, legs cannot attempt establishment yet";
                                        if (!(cause.size() == 3)) throw
                                                "keys in map 'attempt' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        if (cause[0].IsNull()) { // outgoing leg
                                            rat12 = label2rat.at(cause[1].as<string>());
                                            et2 = label2et.at(cause[2].as<string>());
                                            rat23 = NO_RAT;
                                        } else { // incoming leg
                                            rat12 = NO_RAT;
                                            et2 = label2et.at(cause[0].as<string>());
                                            rat23 = label2rat.at(cause[1].as<string>());
                                            if (!(cause[2].IsNull())) throw
                                                    "keys in map 'attempt' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        }
                                    }
                                    inflt2attempt_rate[{ evt, { rat12, et2, rat23 } }] = ar;
                                    inflt2delta_probunits[{ evt, { rat12, et2, rat23 } }] = 0.0;
                                } else { // basic
                                    if ((cause.as<string>() != "basic") && (cause.as<string>() != "base")) throw
                                            "keys in map 'attempt' can be 'basic', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                    evt2base_attempt_rate[evt] = ar;
                                }
                            }
                        }
                    }
                    n3 = n2["success"];
                    if (n3) {
                        if (!n3.IsMap()) {
                            if (!n3.IsScalar()) throw "yaml field 'success' within 'dynamics' must be a scalar (base probunits) or a map";
                            evt2base_probunits[{ ec, et1, rat13, et3 }] = parse_double(n3.as<string>());
                        } else {
                            for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                                auto cause = it3->first;
                                entity_type et2;
                                relationship_or_action_type rat12, rat23;
                                if (cause.IsSequence()) {
                                    auto pu = parse_double(it3->second.as<string>());
                                    if (cause.size() == 5) { // angle
                                        if (!((cause[0].IsNull() || (cause[0].as<string>() == et1l)) && (cause[4].IsNull() || (cause[4].as<string>() == et3l)))) throw
                                                "keys in map 'success' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        rat12 = label2rat.at(cause[1].as<string>());
                                        et2 = label2et.at(cause[2].as<string>());
                                        rat23 = label2rat.at(cause[3].as<string>());
                                    } else {
                                        // TODO later
                                        throw "sorry, legs cannot influence establishment success yet";
                                        if (!(cause.size() == 3)) throw
                                                "keys in map 'success' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        if (cause[0].IsNull()) { // outgoing leg
                                            rat12 = label2rat.at(cause[1].as<string>());
                                            et2 = label2et.at(cause[2].as<string>());
                                            rat23 = NO_RAT;
                                        } else { // incoming leg
                                            rat12 = NO_RAT;
                                            et2 = label2et.at(cause[0].as<string>());
                                            rat23 = label2rat.at(cause[1].as<string>());
                                            if (!(cause[2].IsNull())) throw
                                                    "keys in map 'success' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        }
                                    }
                                    inflt2delta_probunits[{ evt, { rat12, et2, rat23 } }] = pu;
                                } else if ((cause.as<string>() == "basic") || (cause.as<string>() == "base")) {
                                    evt2base_probunits[{ ec, et1, rat13, et3 }] = parse_double(it3->second.as<string>());
                                } else if (cause.as<string>() == "tails") {
                                    if (it3->second.IsSequence()) { // separate tail indices
                                        if (!(it3->second.size() ==2 )) throw
                                                "tail specification must be either a number or a pair of numbers.";
                                        evt2left_tail[evt] = parse_double(it3->second[0].as<string>());
                                        evt2right_tail[evt] = parse_double(it3->second[1].as<string>());
                                    } else { // same tail index
                                        evt2left_tail[evt] = evt2right_tail[evt] = parse_double(it3->second.as<string>());
                                    }
                                    if ((!(evt2left_tail[evt]>=0.0)) || (!(evt2right_tail[evt]>=0.0)) || (!(evt2left_tail[evt]<INFINITY)) || (!(evt2right_tail[evt]<INFINITY))) throw
                                            "tail indices must be non-negative finite numbers";
                                } else {
                                    throw "keys in map 'success' can be 'tails', 'basic', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                }
                            }
                        }
                    }
                }
                n2 = spec["terminate"];
                ec = EC_TERM;
                if (n2) {
                    event_type evt = { ec, et1, rat13, et3 };
                    if (!n2.IsMap()) throw "yaml field 'terminate' within 'dynamics' must be a map";
                    auto n3 = n2["attempt"];
                    if (n3) {
                        evt2left_tail[evt] = evt2right_tail[evt] = TAIL_DEFAULT;
                        evt2base_probunits[evt] = 0.0;
                        if (!n3.IsMap()) {
                            if (!n3.IsScalar()) throw "yaml field 'attempt' within 'dynamics' must be a scalar (base attempt rate) or a map";
                            evt2base_attempt_rate[evt] = parse_double(n3.as<string>());
                        } else {
                            for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                                auto cause = it3->first;
                                auto ar = parse_double(it3->second.as<string>());
                                if (!(ar >= 0.0)) throw "values in map 'attempt' must be non-negative";
                                entity_type et2;
                                relationship_or_action_type rat12, rat23;
                                if (cause.IsSequence()) {
                                    if (cause.size() == 5) { // angle
                                        if (!((cause[0].IsNull() || (cause[0].as<string>() == et1l)) && (cause[4].IsNull() || (cause[4].as<string>() == et3l)))) throw
                                                "keys in map 'attempt' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        rat12 = label2rat.at(cause[1].as<string>());
                                        et2 = label2et.at(cause[2].as<string>());
                                        rat23 = label2rat.at(cause[3].as<string>());
                                    } else {
                                        // TODO later
                                        throw "sorry, legs cannot attempt termination yet";
                                        if (!(cause.size() == 3)) throw
                                                "keys in map 'attempt' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        if (cause[0].IsNull()) { // outgoing leg
                                            rat12 = label2rat.at(cause[1].as<string>());
                                            et2 = label2et.at(cause[2].as<string>());
                                            rat23 = NO_RAT;
                                        } else { // incoming leg
                                            rat12 = NO_RAT;
                                            et2 = label2et.at(cause[0].as<string>());
                                            rat23 = label2rat.at(cause[1].as<string>());
                                            if (!(cause[2].IsNull())) throw
                                                    "keys in map 'attempt' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        }
                                    }
                                    inflt2attempt_rate[{ evt, { rat12, et2, rat23 } }] = ar;
                                    inflt2delta_probunits[{ evt, { rat12, et2, rat23 } }] = 0.0;
                                } else { // basic
                                    if ((cause.as<string>() != "basic") && (cause.as<string>() != "base")) throw
                                            "keys in map 'attempt' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                    evt2base_attempt_rate[evt] = ar;
                                }
                            }
                        }
                    }
                    n3 = n2["success"];
                    if (n3) {
                        if (!n3.IsMap()) {
                            if (!n3.IsScalar()) throw "yaml field 'success' within 'dynamics' must be a scalar (base probunits) or a map";
                            evt2base_probunits[{ ec, et1, rat13, et3 }] = parse_double(n3.as<string>());
                        } else {
                            for (YAML::const_iterator it3 = n3.begin(); it3 != n3.end(); ++it3) {
                                auto cause = it3->first;
                                entity_type et2;
                                relationship_or_action_type rat12, rat23;
                                if (cause.IsSequence()) {
                                    auto pu = parse_double(it3->second.as<string>());
                                    if (cause.size() == 5) { // angle
                                        if (!((cause[0].IsNull() || (cause[0].as<string>() == et1l)) && (cause[4].IsNull() || (cause[4].as<string>() == et3l)))) throw
                                                "keys in map 'success' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        rat12 = label2rat.at(cause[1].as<string>());
                                        et2 = label2et.at(cause[2].as<string>());
                                        rat23 = label2rat.at(cause[3].as<string>());
                                    } else {
                                        // TODO later
                                        throw "sorry, legs cannot influence termination success yet";
                                        if (!(cause.size() == 3)) throw
                                                "keys in map 'success' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        if (cause[0].IsNull()) { // outgoing leg
                                            rat12 = label2rat.at(cause[1].as<string>());
                                            et2 = label2et.at(cause[2].as<string>());
                                            rat23 = NO_RAT;
                                        } else { // incoming leg
                                            rat12 = NO_RAT;
                                            et2 = label2et.at(cause[0].as<string>());
                                            rat23 = label2rat.at(cause[1].as<string>());
                                            if (!(cause[2].IsNull())) throw
                                                    "keys in map 'success' can be 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                        }
                                    }
                                    inflt2delta_probunits[{ evt, { rat12, et2, rat23 } }] = pu;
                                } else if ((cause.as<string>() == "basic") || (cause.as<string>() == "base")) {
                                    evt2base_probunits[evt] = parse_double(it3->second.as<string>());
                                } else if (cause.as<string>() == "tails") {
                                    if (it3->second.IsSequence()) { // separate tail indices
                                        if (!(it3->second.size() ==2 )) throw
                                                "tail specification must be either a number or a pair of numbers.";
                                        evt2left_tail[evt] = parse_double(it3->second[0].as<string>());
                                        evt2right_tail[evt] = parse_double(it3->second[1].as<string>());
                                    } else { // same tail index
                                        evt2left_tail[evt] = evt2right_tail[evt] = parse_double(it3->second.as<string>());
                                    }
                                    if ((!(evt2left_tail[evt]>=0.0)) || (!(evt2right_tail[evt]>=0.0)) || (!(evt2left_tail[evt]<INFINITY)) || (!(evt2right_tail[evt]<INFINITY))) throw
                                            "tail indices must be positive finite numbers";
                                } else {
                                    throw "keys in map 'success' can be 'tails', 'basic'/'base', [~, rel./act.type, ent.type, rel./act.type, ~], [~, rel./act.type, ent.type], or [ent.type, rel./act.type, ~]";
                                }
                            }
                        }
                    }
                }
            }
        } catch (const std::exception&) {
            throw "some entity or the relationship or action type was not declared";
        }
    }

    n1 = c["visualization"];
    if (n1) {
        if (!n1.IsMap()) throw "yaml field 'visualization' must be a map with key entity type";
        for (YAML::const_iterator it1 = n1.begin(); it1 != n1.end(); ++it1) {
            auto label = (it1->first).as<string>();
            auto n2 = it1->second;
            if (!n2.IsSequence()) throw
                    "yaml value fields under 'visualization' must be type: [size/thickness, shape, red, green, blue, alpha]";
            if (label2et.count(label) > 0) {
                auto et = label2et.at(label);
                et2gexf_size[et] = parse_double(n2[0].as<string>());
                et2gexf_shape[et] = n2[1].as<string>();
                et2gexf_r[et] = parse_int(n2[2].as<string>());
                et2gexf_g[et] = parse_int(n2[3].as<string>());
                et2gexf_b[et] = parse_int(n2[4].as<string>());
                et2gexf_a[et] = parse_double(n2[5].as<string>());
            } else {
                if (label2rat.count(label) == 0) throw "relationship/action type has not been declared";
                auto rat = label2rat.at(label);
                rat2gexf_thickness[rat] = parse_double(n2[0].as<string>());
                rat2gexf_shape[rat] = n2[1].as<string>();
                rat2gexf_r[rat] = parse_int(n2[2].as<string>());
                rat2gexf_g[rat] = parse_int(n2[3].as<string>());
                rat2gexf_b[rat] = parse_int(n2[4].as<string>());
                rat2gexf_a[rat] = parse_double(n2[5].as<string>());
            }
        }
    }

    if (!quiet) cout << "...READING CONFIG FINISHED." << endl << endl;
}

