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
* optionally run the tests: ``python3 tests/run_tests.py`` (needs only the Python standard library) and the throughput
  benchmark ``python3 tests/benchmark.py``; see ``CONTRIBUTING.md`` for the repository layout and conventions

Usage
-----
* write some config file ``someconfigfile.yaml`` (see below)
* run model with ``tricl someconfigfile.yaml [options]`` (or first list options with ``tricl someconfigfile.yaml --help``)
* make a movie of the resulting temporal network with ``python/tricl_movie.py`` (see below), or visualize or analyse
  the output gexf-file, e.g. with gephi <https://gephi.org/>

Command line options:
* ``--seed N``: random seed (0 = choose randomly), overrides ``options:seed``
* ``--quiet``, ``--verbose``, ``--debug``, ``--silent``: amount of output, override the corresponding config file options
* ``--logl``: output only the final log-likelihood
* ``--summary``: output only a one-line JSON summary of the final state (number of events, model time, log-likelihood, numbers of links and angles, total event rate, seed used)
* ``--events-out FILE``: write all performed events to a csv file with columns ``t,event,source,relationship,target`` (overrides ``files:events``)
* ``--events-in FILE``: replay the events in the csv file instead of simulating, and compute their log-likelihood (see below)
* ``--grad``: also compute the gradient of the log-likelihood with respect to all model parameters (output with ``--summary``)
* ``--dump-parameters``: only output the model parameters and their current values as JSON and exit
* ``--dump-model``: only output the model structure (entity types with counts, relationship types, link types with initial link counts, event types with influences) as JSON after initialization and exit
* ``--stats-out FILE`` and ``--stats-every DT``: write the numbers of links by link type (and the total numbers of links and angles) to a csv file every ``DT`` model time units (also ``files:stats`` in the config file)
* ``--links-out FILE``: write the time interval of every link to a csv file with columns ``source,relationship,target,start,end`` (the temporal network, see below; also ``files:links``)
* ``--entities-out FILE``: write all entities to a csv file with columns ``id,label,type`` (also ``files:entities``)
* ``--max-events N``, ``--max-t T``: override ``limits:events`` and ``limits:t``; ``--max-wall SECONDS`` (also ``limits:wall``): stop the simulation after this much wall-clock time (the run then ends with the state at the last event, like with the event limit)
* ``--NAME VALUE`` (or ``-X VALUE`` for one-letter names): override the metaparameter ``NAME`` defined in the config file by a value or expression

Caution: output files might get large! Try with small ``limits:events`` first and use gexf.gz file format!

Temporal network output and movies
----------------------------------
``tricl config.yaml --links-out links.csv --entities-out entities.csv`` writes the simulated temporal network as an
*interval list*: one row ``source,relationship,target,start,end`` per directed link and time interval during which it
existed (rows are written when a link is terminated, and at the end of the run for all links still existing, with
the final model time as ``end``), plus one row ``id,label,type`` per entity. Symmetric relationship types give two
rows per pair of entities, one per direction (filter ``source < target`` for an undirected representation);
relationship types that were only declared as the inverse of another type are not written since they are implied.
This plain table is the usual input form of temporal network libraries that read edge lists with times (e.g.
pathpy, teneto, DyNetx, Raphtory, typically via pandas). It is also written in replay mode, so a movie can be made of
an observed event sequence.

``python/tricl_movie.py entities.csv links.csv --out movie.mp4`` renders the temporal network as a movie
(mp4/webm/mkv via ``ffmpeg``, animated gif via Pillow, or a directory of png frames; needs numpy and matplotlib).
Node positions come from a force-directed layout of the *time-aggregated* network (pairs are attracted in proportion
to the total time they were linked), so that the picture stays still and only the links and states move;
``--layout dynamic`` recomputes the layout for every frame instead. Two links in opposite directions are drawn as
one line, a single directed link as an arrow. Nodes are coloured by entity type, or with ``--state R [R ...]`` by
the entity's current link of relationship type R, which is how tricl models usually encode an entity's state
(a link to a "hub" entity such as ``active`` or ``covid-19``); the hub entities are then not drawn. ``--draw``
selects the relationship types to draw, ``--hide`` the entity types not to draw, ``--t0``/``--t1``/``--frames``/``--dt``
the time window and resolution; see ``--help``. Examples:

    tricl config_files/granovetter_helfmann.yaml --seed 1 --E 5000 --entities-out entities.csv --links-out links.csv
    python3 python/tricl_movie.py entities.csv links.csv --state is "is not" --out granovetter.mp4

    tricl config_files/sir_sd.yaml --seed 1 --entities-out entities.csv --links-out links.csv
    python3 python/tricl_movie.py entities.csv links.csv --draw "regularly meets" \
        --state "is susceptible to" "is infectious for" "has recovered from" "has died of" --out sir.gif

The gexf output and the Gephi-based movie workflow in ``gephi/`` remain available.

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

Macroscopic approximation
-------------------------
``python3 python/tricl_macro.py config.yaml [--seeds 1 2 3] [--t-max T] [--dt DT] [--closure poisson|mean]``
derives a system of ordinary differential equations for the numbers of links of every link type from the model
(read via ``tricl config.yaml --dump-model``), integrates it, and optionally compares it with simulations for the
given seeds (using ``--stats-out``, which writes the numbers of links by link type every ``--stats-every`` model
time units). The closure treats the numbers of angles adjacent to a dyad as independent Poisson variables whose
means are products of the current link densities, and averages the event rates over them; it therefore neglects
all correlations between links, and its deviation from the simulation measures the effect of network structure
(clustering, triadic closure). For a model with constant rates and no angle influences, the approximation is exact
in expectation (see ``tests/configs/mean_field.yaml``). Rules of the form "terminate a link immediately when a link
of another type exists on the same dyad" (the idiom for mutually exclusive states, e.g. ``is``/``is not``) are
recognised and handled as state switches; other immediate events (infinite rates) are replaced by a large finite
rate, so models relying on them are only approximated roughly.

RDF and knowledge graphs
------------------------
The tricl data model maps almost one-to-one onto RDF: entities are resources, entity types are classes, relationship
types are object properties (symmetric ones ``owl:SymmetricProperty``, named inverses ``owl:inverseOf``), links are
triples, and angles are property paths of length two. ``python/tricl_rdf.py`` converts in both directions:

* ``python3 python/tricl_rdf.py rdf2tricl DATA.nt --out DIR`` reads RDF data (N-Triples natively, other formats via
  the optional ``rdflib`` package, or a SPARQL endpoint URL with ``--query``) and writes ``DIR/model.yaml``, a config
  skeleton with all entity types, entities, relationship types and ``initial links`` entries, plus one csv file per
  link type. Literal-valued properties are ignored, self-links dropped, symmetric duplicates and inverse pairs folded.
  Add a ``dynamics`` section and run tricl from within ``DIR``. See ``tricl_rdf.py --help`` for options mapping
  classes and properties to labels.
* ``python3 python/tricl_rdf.py tricl2rdf OUTPUT.gexf.gz --out FILE.ttl`` converts the temporal network written by
  tricl into RDF-star (Turtle-star): every interval during which a link existed becomes
  ``<< :source :relationship :target >> tricl:from T1 ; tricl:until T2 .``

A first fit to real data: military alliances
--------------------------------------------
``python/tricl_atop.py --fit`` estimates a model of the formation and dissolution of defense pacts between states
from the ATOP alliance data (Leeds et al. 2002, v5.1) and the Correlates of War state system membership list, as
shipped in the R package `peacesciencer` (downloaded from its GitHub repository; needs ``pyreadr``). Every maximal
run of years in which a pair of states has a defense pledge becomes one link interval of a symmetric relationship
type ``defense``; membership of a state in the international system is a link ``[state, is in, system]`` to a hub
entity, so that "both states are members" is the angle ``[~, is in, system, contains, ~]``. The model has six
metaparameters: formation attempts at rate ``a0`` per pair of member states plus ``a1`` per common defense ally
(the angle ``[~, defense, state, defense, ~]``), dissolution attempts at rate ``d0`` with success probability units
``-b1`` per common ally (``tails: 0``), and entry and exit rates of states. Data: 217 states, 1816-2018, 3000 pact
intervals (1707 still in effect in 2018), 4514 events. Maximum-likelihood estimates (standard errors from the
observed information; the whole fit takes about a minute):

| parameter | estimate | meaning |
|---|---|---|
| ``a0`` | 1.78e-4 ± 0.08e-4 /yr | formation rate per pair of member states without common allies |
| ``a1`` | 1.02e-2 ± 0.02e-2 /yr | additional formation rate per common defense ally |
| ``d0`` | 0.108 ± 0.005 /yr | dissolution attempt rate |
| ``b1`` | 0.037 ± 0.001 | probability units by which each common ally lowers the dissolution success |
| ``entry``, ``exit`` | 0.0081 ± 0.0005, 0.0031 ± 0.0004 /yr | entry and exit rates of states |

So a pair with one common ally forms a pact about 58 times as fast as a pair without, and the dissolution hazard
``d0 · expit(-4 b1 k)`` falls from 0.054/yr without common allies (a mean lifetime of 18 years) to 0.035/yr with
five, 0.020/yr with ten and 0.002/yr with 26 (the size of NATO). Dropping the ally-of-ally term lowers the
log-likelihood from -20841 to -30155, dropping the stabilisation term to -21618. Two cautions: the likelihood in
``b1`` has a second, much worse local maximum at negative values (a saturated sigmoid), so the fit must start near
the optimum (the generated config does; ``tricl_fit.py`` from other start values can end there); and a multilateral
treaty appears as many simultaneous dyadic events, whose order within a year is a convention (by state codes;
with a random order, ``--shuffle``, ``a0`` changes by a quarter and ``a1`` and ``b1`` by a few percent), so both
effects partly reflect treaty-level events. Other pledge types (``--types``), contiguity, disputes or major-power
status can be added as further relationship types and hub entities in the same way.

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
    stats: <csv file to write the numbers of links by link type to at regular model time intervals>
    links: <csv file to write the time interval of every link to>  # columns: source, relationship, target, start, end
    entities: <csv file to write all entities to>  # columns: id, label, type
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
    wall: <max. wall-clock time of the simulation in seconds>  # optional
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
                tails: <common tail index>  # or [<left tail index>, <right tail index>], default: 1.0 (see below)
                base: <basic success probability units>  # default: 0.0
                [~, <relationship label>, <entity type label>, <relationship label>, ~]: <additional success probability unist due to this type of angle>  # may be negative
        terminate:  # specify how this type of link is terminated (if at all)
            # specify attempt, success as above
```
Success probability units ``pu`` are converted into a success probability by a sigmoidal function with two tail indices
``σ0`` (left tail) and ``σ1`` (right tail):
```
f(pu) = T_σ0(-v) / 2 + 1/2 - T_σ1(v) / 2,   v = 4 pu / (k(σ0) + k(σ1)),
T_σ(x) = (1 + σ ln(1 + e^x))^(-1/σ) for σ > 0,   T_0(x) = 1 / (1 + e^x),
k(σ) = (1 + σ ln 2)^(-1 - 1/σ),   k(0) = 1/2.
```
For tail indices ``[0, 0]`` this is exactly the expit (logistic) function of ``4 pu``, ``1 / (1 + e^-4pu)``. A positive tail
index makes the corresponding tail decay like a power law with exponent ``-1/σ`` instead of exponentially (the tails of
the expit function are replaced by the q-exponential ``(1 + σy)^(-1/σ)`` of ``y = ln(1 + e^x)``). For all tail indices the
function is strictly increasing, depends continuously on the tail indices, and has slope 1 at ``pu = 0`` (``k(σ)`` is
minus twice the slope of ``T_σ`` at 0), so near zero, probability units are probability differences. The default tail
index of 1 gives tails ``~ 1/|pu|``, and for tail indices ``[1, 1]`` the function is the same as in earlier versions of tricl.

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
- new form of the sigmoidal function converting probability units to probabilities (q-exponential tails, see above):
  the old formula was broken for exactly one zero tail index (that tail was frozen at 1, so success probabilities
  never fell below 1/2) and its scale constant for tail index 0 was 0 due to integer division; the new one is
  continuous in the tail indices and normalised to slope 1 at zero probability units for *all* tail indices, as the
  old general formula was. For the default tail indices ``[1, 1]`` it is the same function as before, so configs
  using the default tails are unaffected. For tail indices ``[0, 0]`` it is the expit function of ``4 pu`` (the old
  code used the expit of ``pu`` there, with slope 1/4), so configs with ``tails: 0`` and finite success probability
  units (``granovetter_simple.yaml``, ``granovetter_helfmann.yaml``) now have a 4 times steeper sigmoid; their
  probability units of ±10 still give probabilities of practically 0 and 1 as intended, so their behaviour changes
  only marginally (divide the units by 4 to restore it exactly). The old formula also overflowed to a probability
  of exactly 0 below about -124 probability units with the default tails, where the power-law tail actually gives
  small positive values (e.g. 8.7e-5 at -1000, as used in ``sir_sd.yaml``, whose results therefore change; use ``-inf``
  or ``tails: 0`` for events that should be practically impossible). Unit tests in ``tests/test_sigmoid.cpp``.
- fixed: random initial links of a symmetric relationship type between entities of *different* types were never
  generated when the entity ids of the source type happened to be larger than those of the target type (which
  depended on the internal ordering of the entity types). In ``granovetter_helfmann.yaml`` this meant that there were
  no ``knows`` links between always active, contingent and never active agents at all; results of that config change.
- ``python/tricl_atop.py``: a first fit to real data (formation and dissolution of military alliances, see above)
- new options ``--max-events``, ``--max-t`` and ``--max-wall`` (also ``limits:wall``); ``CONTRIBUTING.md``, a throughput
  benchmark ``tests/benchmark.py`` and a sanitizer build in the continuous integration
- new options ``--links-out`` and ``--entities-out`` (also ``files:links``, ``files:entities``) writing the temporal
  network as a plain interval list, and ``python/tricl_movie.py`` rendering it as a movie (see above); the Gephi
  workflow in ``gephi/`` is now legacy
- ``python/tricl_macro.py`` generates and integrates a mean-field approximation of a model and compares it with
  simulations; new options ``--dump-model``, ``--stats-out``, ``--stats-every``
- ``python/tricl_rdf.py`` converts RDF data into config skeletons and gexf output into RDF-star
- fixed: csv quoting of entity labels was kept in the labels; labels were not escaped in gexf output;
  gzipped gexf output was written through a dangling stream pointer (worked by accident) and is now flushed properly
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

