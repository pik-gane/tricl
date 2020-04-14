# tricl
TriCl model in C++

Dependencies
------------
* C++ language standard >=17
* rapidcsv: <https://github.com/d99kris/rapidcsv>

Building
--------
After cloning to some place, e.g. to ``git/tricl``,

```shell
   cd git/tricl
   mkdir -p build/default
   cd build/default
   cmake ../../
   cmake --build .
   cp src/tricl wherever_you_want_the_binary
```

Profiling
---------   
```shell
   valgrind --tool=callgrind --callgrind-out-file=/tmp/callgrind.out tricl; kcachegrind /tmp/callgrind.out
```

Usage
-----
* write some config file ``someconfigfile.yaml`` (see below)
* run model with ``tricl someconfigfile.yaml``
* visualize or analyse output gexf-file, e.g. with gephi <https://gephi.org/>

Caution: output files might get large! Try with small ``limits:events`` first.

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
.inf  # infinity
```

A tricl config file has this overall structure (where stuff in ``<this kind of brackets>`` is a placeholder):
```yaml
metadata: <model and simulation metadata>

files: <where to get and put stuff>
options: <...>

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
    gexf: <where to output the resulting temporal network>.gexf 
    diagram prefix: <filename prefix for structural diagram output>
    # files not listed are not generated

limits:
    t: <max. simulation time>  # default: .inf
    events: <max. no. of simulated events>  # default: .inf

options:
    quiet: <true or false>  # default: false
    verbose: <true or false>  # default: false
    debug: <true or false>  # default: false
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