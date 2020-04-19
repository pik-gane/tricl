# tricl
TriCl model in C++

Dependencies
------------
* C++ language standard >=17
* rapidcsv: <https://github.com/d99kris/rapidcsv>
* zlib <https://www.zlib.net/> (instead of using the heavier and harder to install libboost_iostreams shared library)
* boost::iostreams header files: <https://www.boost.org/doc/libs/1_72_0/libs/iostreams/doc/index.html> 

Installation
------------
* ``cd`` to some place
* ``git clone https://github.com/mensch72/tricl.git``
* ``cd tricl``
* ``mkdir -p build/default``
* ``cd build/default``
* if necessary, set the environmental variables CC, CXX, CPATH, LIBRARY_PATH, LD_LIBRARY_PATH to point to your C compiler, C++ compiler, static library path, an shared library path
* ``cmake ../../``
* ``cmake --build .``
* ``cp src/tricl`` to wherever you want the binary

Usage
-----
* write some config file ``someconfigfile.yaml`` (see below)
* run model with ``tricl someconfigfile.yaml [options]`` (or first list options with ``tricl someconfigfile.yaml --help``)
* visualize or analyse output gexf-file, e.g. with gephi <https://gephi.org/>

Caution: output files might get large! Try with small ``limits:events`` first and use gexf.gz file format!

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
In addition, numerical values can not only be specified as numeric lieterals but also via simple mathematical expressions such as ``3 * sin(pi/5)^2``.

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
    # files not listed are not generated

options:
    quiet: <true or false>  # default: false
    verbose: <true or false>  # default: false
    debug: <true or false>  # default: false

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
?

Change log
----------

2020-04-18
- automatically generate code documentation and publish it to <https://mensch72.github.io/tricl/html/index.html> following this tutorial: <https://gist.github.com/francesco-romano/351a6ae457860c14ee7e907f2b0fc1a5>

2020-04-17
- metaparameters can now be overwritten on the command line

2020-04-16
- added utility scripts for making videos with gephi (folder ``gephi``)
- added support for gzipped output (.gexz.gz, not on the cluster yet) and separate output of individual relationship types

Development
-----------

Code documentation is here: <https://mensch72.github.io/tricl/html/index.html>
(automatically updated after each commit to the master branch).
 
To profile:
```shell
   valgrind --tool=callgrind --callgrind-out-file=/tmp/callgrind.out tricl; kcachegrind /tmp/callgrind.out
```

To generate a local copy of the documentation:
* execute ``doxygen`` in the top repository folder
* find the documentation in the folder ``doc`` (will be ignored by git)
