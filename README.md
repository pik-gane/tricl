# tricl
TriCl model in C++

Dependencies
------------
* C++ language standard >=17 and CMake >= 3.10
* yaml-cpp: <https://github.com/jbeder/yaml-cpp> (Debian/Ubuntu package ``libyaml-cpp-dev``)
* zlib <https://www.zlib.net/> (package ``zlib1g-dev``; used instead of the heavier and harder to install libboost_iostreams shared library)
* boost header files (package ``libboost-dev``): ``boost::container::flat_set`` and ``boost::iostreams``
* optional: graphviz (``dot``) for the structural diagrams, doxygen for the code documentation
* bundled in ``src/3rdparty``: rapidcsv, tinyexpr, cxxopts, and two source files of boost::iostreams

Installation
------------
* ``git clone https://github.com/pik-gane/tricl.git``
* ``cd tricl``
* if necessary, set the environmental variables CC, CXX, CPATH, LIBRARY_PATH, LD_LIBRARY_PATH to point to your C compiler, C++ compiler, static library path, an shared library path
* ``cmake -S . -B build`` (add ``-DCMAKE_BUILD_TYPE=Debug`` for a build with assertions and debug symbols; the default is an optimized ``Release`` build)
* ``cmake --build build``
* ``cp build/src/tricl`` to wherever you want the binary
* optionally run the tests: ``python3 tests/run_tests.py`` (needs only the Python standard library)

Usage
-----
* write some config file ``someconfigfile.yaml`` (see below)
* run model with ``tricl someconfigfile.yaml [options]`` (or first list options with ``tricl someconfigfile.yaml --help``)
* visualize or analyse output gexf-file, e.g. with gephi <https://gephi.org/>

Command line options:
* ``--seed N``: random seed (0 = choose randomly), overrides ``options:seed``
* ``--quiet``, ``--verbose``, ``--debug``, ``--silent``: amount of output, override the corresponding config file options
* ``--logl``: output only the final log-likelihood
* ``--summary``: output only a one-line JSON summary of the final state (number of events, model time, log-likelihood, numbers of links and angles, total event rate, seed used)
* ``--events-out FILE``: write all performed events to a csv file with columns ``t,event,source,relationship,target`` (overrides ``files:events``)
* ``--events-in FILE``: replay the events in the csv file instead of simulating, and compute their log-likelihood (see below)
* ``--grad``: also compute the gradient of the log-likelihood with respect to all model parameters (output with ``--summary``)
* ``--dump-parameters``: only output the model parameters and their current values as JSON and exit
* ``--NAME VALUE`` (or ``-X VALUE`` for one-letter names): override the metaparameter ``NAME`` defined in the config file by a value or expression

Caution: output files might get large! Try with small ``limits:events`` first and use gexf.gz file format!

Log-likelihood, replay mode and parameter estimation
----------------------------------------------------
Each run reports the log-likelihood of the simulated trajectory (``logl``): the sum over events of the log of the
event's rate, minus the integral of the total event rate over time (the log-probability that nothing else happened),
including the final interval until ``limits:t`` if that is finite.

The same quantity can be computed for a *given* sequence of events: ``tricl config.yaml --events-in events.csv --summary``
replays the events in the csv file (columns ``t,event,source,relationship,target``, as written by ``--events-out``;
``event`` is ``establish`` or ``terminate``) instead of simulating, starting from the initial links in the config file,
and reports their log-likelihood under the model. The sequence must be complete, i.e. it must also contain the
"immediate" events (with infinite rates) that the model performs in response to other events; events that are
impossible under the model are reported as errors.

Replaying draws no random numbers and keeps no schedule. The observation window ends at ``limits:t`` (the
log-probability of no further event until then is included), or with the last event if ``limits:events`` is reached,
exactly as in a simulation. Random initial links (block and geometric models) are generated from a separate random
number generator seeded from ``--seed``, so the initial state of a simulation and of a replay with the same config
and seed are identical, and replaying the events of a simulation reproduces its log-likelihood.

With ``--grad``, the gradient of the log-likelihood with respect to all model parameters specified in the
``dynamics`` section (base attempt rates, attempt rates and probunits of individual influences) is computed
analytically and included in the ``--summary`` output, keyed by labels like
``"establish that user knows user | attempt via knows user knows"``. ``--dump-parameters`` prints these labels
and their current values without running anything, which allows computing derivatives with respect to
metaparameters by the chain rule. (Parameters can only be estimated if their current value is nonzero, since
angles without effect are not tracked.)

``python/tricl_fit.py config.yaml events.csv --fit NAME ... [--se]`` estimates metaparameters by maximum likelihood
using these facilities (scipy is used if available), optionally with standard errors from the observed information,
and ``python/tricl_recovery.py`` runs a parameter-recovery study (simulate with known values, re-estimate).
See ``tests/configs/three_entities.yaml`` for a config in which every rate is a metaparameter.

Legend to output
----------------
- logl: log-likelihood of this realization so far (log of the probability density of the simulated trajectory given the initial state, including the probability that no other event happened in between; events that happen "immediately" contribute the log-probability of their random order)
- er: current rate of events [1/time]
- ld: overall link density
- ad: overall angle density
- q: rough metric of clustering (= ad/ld²). 1 indicates a "normal level" of clustering
- t: model time
 
Config file syntax
------------------
See folder ``config_files`` for examples.
Config files use YAML <https://yaml.org/> syntax, in particular using these YAML features:
```yaml
simple key: scalar value
simple key: [list, of, values] 
[list, as, key]: value

list name:
   - list item
   - list item

map name:
   key: value
   key: value

# comment
~  # empty value
inf  # infinity
```
In addition, numerical values can not only be specified as numeric literals but also via simple mathematical expressions such as ``3 * sin(pi/5)^2`` (using the syntax of tinyexpr <https://github.com/codeplea/tinyexpr>, extended by the symbols ``inf``, ``infinity``, ``eps``, and the metaparameters defined in the config file). Infinity can be written as ``inf`` or as YAML's ``.inf``. Expressions that cannot be parsed (e.g. because they use an undefined metaparameter) are reported as errors.

A tricl config file has this overall structure (where stuff in ``<this kind of brackets>`` is a placeholder):
```yaml
metadata: <model and simulation metadata>

files: <where to get and put stuff>
options: <...>

metaparameters: <optional definition of values to be used further down>

limits: <how long to run etc.>

entities: <specification of entity types and named entities>
relationship types: <spec. of rel.types and their inverses>

initial links:
    named: <named initial links>
    random: <spec. for random generation of initial links>
    <filename.csv>: <spec. for reading initial links from a file>

dynamics: <rules of the model dynamics>

visualization: <node and link attributes to be used by gephi>  # optional
```
These sections look like this in more detail:
```yaml
metadata:
    name: <name of the model> 
    description: 
        <some lengthier description,
        which might be several lines>
    author: <who wrote this file>
    date: <...>
    email: <...>

files:   
    gexf: <where to output the resulting temporal network>  
        # must end in either .gexf or .gexf.gz (recommended) 
    diagram prefix: <filename prefix for structural diagram output>
    events: <csv file to write all performed events to>  # columns: t, event, source, relationship, target
    # files not listed are not generated

options:
    silent:  <true or false>  # suppress all output, default: false
    logl:    <true or false>  # only output final log-likelihood, default: false
    quiet:   <true or false>  # suppress most output, default: false
    verbose: <true or false>  # increase output, default: false
    debug:   <true or false>  # output very much, default: false

metaparameters:  
    # will be substituted for their values 
    # in expressions occurring in the rest of the config file.
    # a metaparameter XY can be overwritten by specifying 
    # an option -XY <value or expression> on the command line 

    <token>: <value or expression>   # <explanation>
    <token>: <value or expression, maybe using earlier tokens>
    <token>: <value or expression, maybe using earlier tokens>

limits:
    t: <max. simulation time>  # default: .inf
    events: <max. no. of simulated events>  # default: .inf
    # at least one of the two must be finite
```
```yaml
entities:  
    # map of entity type specifications.
    # the following formats are possible:

    <entity type label>: <number of entities of this type>

    <entity type label>:  # list of all entities of this type
        - <entity label>
        - <entity label>

    <entity type label>:
        n: <number of entities of this type>
        named:  # labels for some of these n entities
            - <entity label>
            - <entity label>

    # in some cases, single-entity types are needed, 
    # where entity type and entity have the same label, e.g.
    <singular entity label>:
        - <same singular entity label>

relationship types:
    # map of relationship type specifications.
    # the following formats are possible:
    <relationship label>: symmetric  # a symmetric (undirected) relationship type
    <relationship label>: <inverse relationship label>  # a non-symmetric (directed) relationship type with a named inverse
    <relationship label>: ~  # a non-symmetric (directed) relationship type without any named inverse
    <relationship label>: 
        inverse: <~ or inverse relationship label>
        gexf: <~ or filename>  # suppress or redirect output of this type to separate file
```
```yaml
initial links: 
    
    named:  # list of names initial links
        - [<entity label>, <relationship label>, <entity label>]
        - [<entity label>, <relationship label>, <entity label>]

    random:  
        # map of random link specifications, keyed by link type
        # the following formats are possible:

        # Erdös-Renyi model:
        [<entity type label>, <relationship label>, <entity type label>]:
            probability: <link density as floating point number>
    
        # homogeneous block model:
        [<entity type label>, <relationship label>, <entity type label>]:
            blocks: <no. of blocks of equal expected size>
            within: <within-block link density>
            between: <between-block link density>

        # random geometric model:
        [<entity type label>, <relationship label>, <entity type label>]:
            dimension: <no. of spatial dimensions>
            decay: <rate of exponential decay of link probability with distance>

    # specifications for reading initial links from a file.
    # the following formats are possible:

    # single-link-type file:
    <filename.csv>:  # file extension must be ".csv"
        type: [<entity type label>, <relationship label>, <entity type label>]
        cols: [<source column index>, <target column index>]  # default: [0, 1] 
        skip: <no. of rows to skip at beginning, incl. column heads>  # default: 0
        max: <no. of rows to read>  # default: .inf
        delimiter: "<delimiting character>"  # default: ","

    # multi-link-type file:
    <filename.csv>:  # file extension must be ".csv"
        entity types: [<source entity type label>, <target entity type label>]
        cols: [<source column index>, <relationship column index>, <target column index>]  # default: [0, 1, 2] 
        # skip, max, delimiter as above
```
```yaml
dynamics:
    # map of dynamic rules, keyed by link type.

    [<entity type label>, <relationship label>, <entity type label>]:  
        establish:  # specify how this type of link is established (if at all)
            attempt:  # attempt rates
                base: <basic attempt rate>  # default: 0.0
                [~, <relationship label>, <entity type label>, <relationship label>, ~]: <additional attempt rate due to this type of angle>
            success:  # success probability units
                tails: <common tail index>  # or [<left tail index>, <righ tail index>], default: 1.0
                base: <basic success probability units>  # default: 0.0
                [~, <relationship label>, <entity type label>, <relationship label>, ~]: <additional success probability unist due to this type of angle>  # may be negative
        terminate:  # specify how this type of link is terminated (if at all)
            # specify attempt, success as above
```
```yaml    
visualization:  

    # see https://gephi.org/gexf/1.2draft/gexf-12draft-primer.pdf, p.22
    # (gephi still seems to ignore node shape)

    <entity type label>: [<size>, <shape>, <r>,<g>,<b>,<a>]  # size, a: floats. r,g,b: 0...255. shape: disc (default), diamond, square, triangle 

    <relationship label>: [<thickness>, <shape>, <r>,<g>,<b>,<a>]  # thickness: float. shape: solid (default), dotted, dashed, double
```

Third-party code used
---------------------
(in folder ``3rdparty``)
* cxxopts: <https://github.com/jarro2783/cxxopts>
* tinyexpr: <https://github.com/codeplea/tinyexpr>
* gzip.cpp and zlib.cpp from boost::iostreams

License
-------
GNU General Public License v3.0, see file ``LICENSE``.

Change log
----------

2026-09-17
- random initial links now use a separate random number generator (seeded from the seed), so the initial state does
  not depend on the dynamics; trajectories for a given seed therefore differ from earlier versions
- replay mode (``--events-in``) computing the log-likelihood of a given event sequence, analytic gradients of the
  log-likelihood w.r.t. all model parameters (``--grad``), ``--dump-parameters``, and Python scripts for
  maximum-likelihood estimation of metaparameters and parameter-recovery studies (folder ``python``)
- performance: angles are found with a galloping intersection into a reusable buffer, per-type parameters live in
  small dense tables instead of two 512 MB arrays and hash maps, the sigmoid's scale parameter is precomputed,
  leg sets are vectors, and the status line is throttled. Measured (one core, seed 1): peak memory 1029 MB → 10 MB,
  start-up 0.8 s → 0.01 s, ``granovetter_helfmann.yaml`` with k=10: 63k → 158k events/s at N=100,
  14.7k → 147k at N=1000, 1.7k → 76k at N=10000; dense ``parameters_3blocks.yaml``: 163 s → 101 s.
  Simulated trajectories are unchanged for a given seed.
- fixed the log-likelihood computation (it contained a spurious -log(total rate) term, missed the waiting time
  of rejected summary event attempts and the probability of no further event until the time limit, and was infinite
  whenever an "immediate" event was the only pending one)
- fixed the bookkeeping of the total event rate (the summary event share of the inverse link of a symmetric or inverse
  relationship was never removed; the total is now also protected against floating point drift)
- infinite contributions to attempt rates and success probunits are now counted separately, so that removing one of
  several infinite contributions no longer yields NaN
- expressions in config files that cannot be parsed now produce an error instead of silently becoming NaN;
  YAML's ``.inf`` is accepted; unknown entity, entity type and relationship labels are reported by name;
  missing input files are reported by name; ``limits:t`` may be infinite if ``limits:events`` is finite
- new options ``--summary``, ``--events-out``, ``files:events``; one-letter metaparameters can be given as ``--X value``
- the build now uses CMake build types (optimized ``Release`` build by default), ``find_package`` for the dependencies,
  and returns a non-zero exit code on errors
- added a test suite (``tests/run_tests.py``) and GitHub Actions workflows for building, testing and publishing the
  documentation (replacing the defunct Travis CI setup)
- moved the unfinished ``identity.yaml`` to ``config_files/drafts``

2020-04-25
- add computation and output of log-likelihood (not for prescribed trajectories yet) 

2020-04-23
- improved log output
- in csv input, different prefixes for source and target entity labels can be used

2020-04-22
- added command-line options --quiet, --verbose, --seed X overriding the config file options

2020-04-18
- automatically generate code documentation and publish it to <https://mensch72.github.io/tricl/html/index.html> following this tutorial: <https://gist.github.com/francesco-romano/351a6ae457860c14ee7e907f2b0fc1a5>

2020-04-17
- metaparameters can now be overwritten on the command line

2020-04-16
- added utility scripts for making videos with gephi (folder ``gephi``)
- added support for gzipped output (.gexz.gz, not on the cluster yet) and separate output of individual relationship types

Development
-----------

Tests: ``python3 tests/run_tests.py [--bin build/src/tricl]`` runs all example configs with fixed seeds and reduced
event limits, compares the final state (no. of events, model time, log-likelihood, numbers of links and angles) with the
reference values in ``tests/reference.json``, checks the log-likelihood of a tiny model against an exact formula, and
checks the error handling of the config parser. After an intentional change of the simulated trajectories (e.g. a change
in the order of random draws), regenerate the references with ``--update``. The GitHub Actions workflow in
``.github/workflows/ci.yml`` runs the tests for every push in Release and Debug builds.

Consistency checks: ``tricl myconfig.yaml --debug`` verifies the internal data structures after every event (slow).

To generate a local copy of the documentation:
* execute ``doxygen`` in the top repository folder
* find the documentation in the folder ``doc`` (will be ignored by git)
* the workflow in ``.github/workflows/docs.yml`` publishes the documentation to the ``gh-pages`` branch on every push to ``master``

To profile:
```shell
   valgrind --tool=callgrind --callgrind-out-file=callgrind.out tricl myconfig.yaml --quiet && kcachegrind callgrind.out &
```

