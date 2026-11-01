Contributing to tricl
=====================

Thanks for considering a contribution. This page explains how the repository is organised, how to build and test
it, and what a good pull request looks like.

Building and testing
--------------------

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release   # or Debug (assertions and --debug consistency checks)
cmake --build build -j 4
python3 tests/run_tests.py                        # all tests, needs only the Python standard library (a few tests
                                                  # need numpy and matplotlib and are skipped otherwise)
python3 tests/run_tests.py --only sigmoid         # a subset by name
python3 tests/benchmark.py                        # throughput of the example configs, see below
```

Dependencies: a C++17 compiler, CMake ≥ 3.16, yaml-cpp, Boost headers, zlib; graphviz for the structural
diagrams. On Debian/Ubuntu: `apt-get install cmake g++ libyaml-cpp-dev libboost-dev zlib1g-dev graphviz`.

The continuous integration (`.github/workflows/ci.yml`) builds Release and Debug, runs the test suite, runs one
example under the Debug consistency checks, and builds with the address and undefined-behaviour sanitizers.

What the tests check
--------------------

`tests/run_tests.py` runs, among others:

* unit tests of the sigmoidal function (`tests/test_sigmoid.cpp`, built by CMake and also available via `ctest`);
* regression tests: every example config with a fixed seed against `tests/reference.json` (number of events,
  log-likelihood, link and angle counts). The references are specific to the random number generator of
  libstdc++; regenerate them with `--update` after a change that legitimately alters trajectories, and say so
  in the pull request;
* an exact log-likelihood test on a two-entity model, replay consistency (simulate, replay, compare the
  log-likelihoods), analytic gradients against finite differences, a maximum-likelihood recovery, the RDF and
  temporal-network exports, the macroscopic approximation, and the error handling of malformed configs.

A change to the dynamics (rates, likelihood, scheduling) needs a test that would have failed before it.

Repository layout
-----------------

* `src/` the simulator. `data_model.h` defines the types, `global_variables.h` declares the global state,
  `config.cpp` reads the config file and command line, `init.cpp` sets up types, entities, initial links and
  event rates, `event.cpp`/`event.h` schedule and perform events and maintain the rates and their gradients,
  `angle.cpp` finds angles, `probability.h` holds the sigmoidal function, `simulate.cpp` the main loop,
  `replay.cpp` the replay of observed events, `io.cpp` and `gexf.cpp` the outputs, `finish.cpp` the end of a run.
  Third-party single-header libraries live in `src/3rdparty`.
* `config_files/` example models; `config_files/drafts/` unfinished ones.
* `python/` tools that drive the binary: parameter estimation (`tricl_fit.py`, `tricl_recovery.py`), RDF
  conversion (`tricl_rdf.py`), the mean-field approximation (`tricl_macro.py`), movies (`tricl_movie.py`), and
  the alliance-data example (`tricl_atop.py`). They need only the standard library unless documented otherwise.
* `tests/` the test runner, reference values, small test configs and data.
* `gephi/` the legacy Gephi movie workflow.

Conventions
-----------

* C++17, no exceptions for control flow except the `throw "message"` convention for config errors, which the
  main function reports. Per-type data live in dense tables indexed by type ids (see `init_types()`); avoid
  hash lookups in the hot loops of `event.cpp` and `angle.cpp`.
* Every command line option has a `files:`/`options:` counterpart in the config file where that makes sense, and
  is listed in the README's options list.
* Config files are the model description language: new syntax needs a README entry in the "config file
  syntax" section and an example config or test config that uses it.
* Document user-visible changes in the README's change log.
* Keep the log-likelihood exact: any change to rates must keep `tests/run_tests.py --only "exact"` and the
  replay/gradient tests passing.

Adding a model
--------------

Copy `config_files/parameters_TEMPLATE.yaml`, describe entity types, relationship types (with inverses or
`symmetric`), initial links and the `dynamics` rules. States of entities are usually encoded as links to a
single hub entity (e.g. `[agent, is, active]`), which makes state-dependent rates expressible as angle
influences. Run with a small `limits: events` first, then add the config to `REGRESSION_CASES` in
`tests/run_tests.py` and regenerate the references with `--update` if it should be covered by the tests.

Pull requests
-------------

Small, focused pull requests with a test and a change-log line are the easiest to review. If a change alters
simulation results, explain why (bug fix, changed semantics) and update the references in the same pull request.
