#!/usr/bin/env python3
"""Exchange of network data between tricl and RDF (Resource Description Framework).

The tricl data model maps almost one-to-one onto RDF: entities are resources, entity types are classes
(rdf:type), relationship types are object properties, links are triples, symmetric relationship types are
owl:SymmetricProperty, and named inverses are owl:inverseOf. Angles are SPARQL property paths of length two.

Usage:

    python3 python/tricl_rdf.py rdf2tricl INPUT [INPUT ...] --out DIR [options]

        Reads RDF data and writes a tricl config skeleton DIR/model.yaml (entity types with all their entities,
        relationship types with symmetry and inverses, and 'initial links' entries) plus one csv file per link
        type. The 'dynamics' section is left for you to fill in. Run tricl from within DIR (the csv paths are
        relative), or edit the paths.

        INPUT may be a file in N-Triples format (parsed natively), any RDF format readable by rdflib if that
        package is installed (Turtle, RDF/XML, JSON-LD, ...), or a SPARQL endpoint URL together with --query.

        Options:
          --query QUERY_OR_FILE   SPARQL SELECT query returning variables ?s ?p ?o (used with an endpoint URL)
          --format FORMAT         rdflib format name for the input files (default: guessed from the extension)
          --default-type LABEL    entity type for resources without an rdf:type (default: "thing")
          --type-map IRI=LABEL    map a class IRI to an entity type label (repeatable; also chooses the type for
                                  resources with several classes: the first mapped class wins)
          --relationship-map IRI=LABEL   map a property IRI to a relationship type label (repeatable)
          --symmetric LABEL       treat a relationship type as symmetric (repeatable; owl:SymmetricProperty statements
                                  in the data are recognised automatically)
          --inverse LABEL1=LABEL2 declare two relationship types as inverses (repeatable; owl:inverseOf is recognised)
          --language LANG         preferred language of rdfs:label values (default: en)
          --name NAME             model name for the metadata section

    python3 python/tricl_rdf.py tricl2rdf OUTPUT.gexf[.gz] --out FILE.ttl [--base IRI]

        Converts the temporal network written by tricl (gexf, optionally gzipped) into RDF-star in Turtle syntax
        ("Turtle-star"): nodes become resources typed by their entity type, and every interval during which a link
        existed becomes an annotated triple << :source :relationship :target >> tricl:from T1 ; tricl:until T2 .

Only the Python standard library is required; rdflib is optional.
"""

import argparse
import csv
import gzip
import json
import os
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

RDF_TYPE = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
RDFS_LABEL = "http://www.w3.org/2000/01/rdf-schema#label"
OWL_SYMMETRIC = "http://www.w3.org/2002/07/owl#SymmetricProperty"
OWL_INVERSE_OF = "http://www.w3.org/2002/07/owl#inverseOf"
VOCABULARY_PREFIXES = ("http://www.w3.org/1999/02/22-rdf-syntax-ns#", "http://www.w3.org/2000/01/rdf-schema#",
                       "http://www.w3.org/2002/07/owl#")
TRICL_NS = "https://github.com/pik-gane/tricl#"


# ---------------------------------------------------------------------------------------------------------------
# reading RDF

class Literal(str):
    """A literal value (as opposed to an IRI or blank node), with optional language tag."""
    lang = ""


NT_LINE = re.compile(r'^\s*(<[^>]*>|_:\S+)\s+(<[^>]*>)\s+(<[^>]*>|_:\S+|"(?:[^"\\]|\\.)*"(?:@[\w-]+|\^\^<[^>]*>)?)\s*\.\s*$')


def _nt_term(token):
    if token.startswith("<"):
        return token[1:-1]
    if token.startswith("_:"):
        return token
    m = re.match(r'^"((?:[^"\\]|\\.)*)"(?:@([\w-]+))?(?:\^\^<[^>]*>)?$', token)
    value = m.group(1).encode("utf-8").decode("unicode_escape")
    lit = Literal(value)
    lit.lang = m.group(2) or ""
    return lit


def read_ntriples(path):
    """Yield (s, p, o) from an N-Triples file (o is a str IRI/blank node or a Literal)."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            m = NT_LINE.match(line)
            if not m:
                raise ValueError("cannot parse line %d of %s: %s" % (n, path, line.strip()[:100]))
            yield _nt_term(m.group(1)), _nt_term(m.group(2)), _nt_term(m.group(3))


def read_rdflib(path, fmt=None):
    """Yield (s, p, o) from any RDF file that rdflib can parse."""
    try:
        import rdflib
    except ImportError:
        raise SystemExit("reading %s needs the rdflib package (pip install rdflib); N-Triples files (.nt) are parsed natively" % path)
    g = rdflib.Graph()
    g.parse(path, format=fmt)
    for s, p, o in g:
        if isinstance(o, rdflib.Literal):
            lit = Literal(str(o))
            lit.lang = o.language or ""
            o = lit
        else:
            o = str(o)
        yield str(s), str(p), o


def read_sparql(endpoint, query):
    """Yield (s, p, o) from the bindings of a SPARQL SELECT query with variables ?s ?p ?o."""
    if os.path.exists(query):
        with open(query) as f:
            query = f.read()
    url = endpoint + ("&" if "?" in endpoint else "?") + urllib.parse.urlencode({"query": query})
    req = urllib.request.Request(url, headers={"Accept": "application/sparql-results+json", "User-Agent": "tricl_rdf.py"})
    with urllib.request.urlopen(req) as resp:
        data = json.load(resp)
    for b in data["results"]["bindings"]:
        if not all(v in b for v in ("s", "p", "o")):
            raise ValueError("the query must return the variables ?s ?p ?o")
        o = b["o"]
        if o["type"] == "literal" or o["type"] == "typed-literal":
            lit = Literal(o["value"])
            lit.lang = o.get("xml:lang", "")
            ov = lit
        else:
            ov = o["value"]
        yield b["s"]["value"], b["p"]["value"], ov


def read_triples(inputs, query=None, fmt=None):
    for inp in inputs:
        if inp.startswith("http://") or inp.startswith("https://"):
            if not query:
                raise SystemExit("a SPARQL endpoint URL needs --query")
            yield from read_sparql(inp, query)
        elif fmt in (None, "nt", "ntriples") and (inp.endswith(".nt") or inp.endswith(".nt.gz")):
            yield from read_ntriples(inp)
        else:
            yield from read_rdflib(inp, fmt)


# ---------------------------------------------------------------------------------------------------------------
# rdf2tricl

def local_name(iri):
    """The part of an IRI after the last '#' or '/', un-percent-encoded."""
    name = re.split(r"[#/]", iri.rstrip("/#"))[-1]
    return urllib.parse.unquote(name) or iri


def yaml_quote(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def unique_labels(items, get_label):
    """Assign unique labels to items (a dict item -> label), appending numbers on clashes."""
    used = {}
    result = {}
    for item in sorted(items):
        label = get_label(item)
        if label in used:
            k = 2
            while "%s (%d)" % (label, k) in used:
                k += 1
            label = "%s (%d)" % (label, k)
        used[label] = item
        result[item] = label
    return result


def rdf2tricl(args):
    triples = []
    types = {}          # resource -> set of class IRIs
    labels = {}         # resource -> {lang: label}
    symmetric = set()   # property IRIs declared owl:SymmetricProperty
    inverse_of = {}     # property IRI -> inverse property IRI
    primary = set()     # of each inverse pair, the property that is kept (the subject of the owl:inverseOf statement)
    n_literal = 0
    for s, p, o in read_triples(args.inputs, args.query, args.format):
        if p == RDF_TYPE and not isinstance(o, Literal):
            if o == OWL_SYMMETRIC:
                symmetric.add(s)
            else:
                types.setdefault(s, set()).add(o)
        elif p == RDFS_LABEL and isinstance(o, Literal):
            labels.setdefault(s, {})[o.lang] = str(o)
        elif p == OWL_INVERSE_OF and not isinstance(o, Literal):
            inverse_of[s] = o
            inverse_of[o] = s
            if o not in primary:
                primary.add(s)
        elif isinstance(o, Literal):
            n_literal += 1
        elif any(p.startswith(v) for v in VOCABULARY_PREFIXES):
            continue
        else:
            triples.append((s, p, o))

    # relationship types:
    rel_map = dict(item.split("=", 1) for item in args.relationship_map)
    properties = sorted(set(p for s, p, o in triples))
    rel_label = unique_labels(properties, lambda p: rel_map.get(p) or local_name(p))
    for label in args.symmetric:
        for p, l in rel_label.items():
            if l == label:
                symmetric.add(p)
    for item in args.inverse:
        l1, l2 = item.split("=", 1)
        p1 = [p for p, l in rel_label.items() if l == l1]
        p2 = [p for p, l in rel_label.items() if l == l2]
        if p1 and p2:
            inverse_of[p1[0]] = p2[0]
            inverse_of[p2[0]] = p1[0]
            primary.add(p1[0])
            primary.discard(p2[0])
        elif p1:
            inverse_of[p1[0]] = "tricl:inverse:" + l2  # inverse without data of its own
            rel_label["tricl:inverse:" + l2] = l2
            primary.add(p1[0])

    # entity types:
    type_map = {}
    for item in args.type_map:
        iri, label = item.split("=", 1)
        type_map[iri] = label
    entities = sorted(set([s for s, p, o in triples] + [o for s, p, o in triples]))
    class_of = {}
    for e in entities:
        classes = types.get(e, set())
        chosen = None
        for iri in type_map:  # first mapped class wins
            if iri in classes:
                chosen = iri
                break
        if chosen is None and classes:
            chosen = sorted(classes)[0]
        class_of[e] = chosen
    class_labels = unique_labels(set(c for c in class_of.values() if c), lambda c: type_map.get(c) or local_name(c))
    et_of = {e: (class_labels[c] if c else args.default_type) for e, c in class_of.items()}

    def best_label(e):
        d = labels.get(e, {})
        for lang in (args.language, "", None):
            if lang in d:
                return d[lang]
        return next(iter(d.values())) if d else local_name(e)

    e_label = unique_labels(entities, best_label)

    # links by link type, with symmetric pairs and inverse properties deduplicated:
    links = {}
    n_self = 0
    seen = set()
    for s, p, o in triples:
        if s == o:
            n_self += 1
            continue
        if p in inverse_of and inverse_of[p] in rel_label and p not in primary:
            # keep only the primary property of an inverse pair, expressed from its side:
            s, p, o = o, inverse_of[p], s
        if p in symmetric:
            key = (p, min(s, o), max(s, o))
        else:
            key = (p, s, o)
        if key in seen:
            continue
        seen.add(key)
        links.setdefault((et_of[s], rel_label[p], et_of[o]), []).append((e_label[s], e_label[o]))

    # write outputs:
    os.makedirs(args.out, exist_ok=True)
    et_entities = {}
    for e, et in et_of.items():
        et_entities.setdefault(et, []).append(e_label[e])
    name = args.name or os.path.splitext(os.path.basename(args.inputs[0]))[0]
    lines = ["metadata:", "    name: " + yaml_quote(name),
             "    description: " + yaml_quote("generated by tricl_rdf.py from " + ", ".join(args.inputs)), "",
             "limits:", "    t: 1  # TODO: simulation time", "",
             "entities:"]
    for et in sorted(et_entities):
        lines.append("    " + yaml_quote(et) + ":")
        for label in sorted(et_entities[et]):
            lines.append("        - " + yaml_quote(label))
    lines += ["", "relationship types:"]
    written_rel = set()
    for p in sorted(rel_label, key=lambda p: (p in inverse_of and p not in primary, rel_label[p])):
        label = rel_label[p]
        if label in written_rel:
            continue
        if p in symmetric:
            lines.append("    " + yaml_quote(label) + ": symmetric")
        elif p in inverse_of and inverse_of[p] in rel_label:
            inv = rel_label[inverse_of[p]]
            lines.append("    " + yaml_quote(label) + ": " + yaml_quote(inv))
            written_rel.add(inv)
        else:
            lines.append("    " + yaml_quote(label) + ": ~")
        written_rel.add(label)
    lines += ["", "initial links:"]
    for (et1, rel, et3), pairs in sorted(links.items()):
        fname = "links_%s_%s_%s.csv" % tuple(re.sub(r"[^A-Za-z0-9]+", "_", x) for x in (et1, rel, et3))
        with open(os.path.join(args.out, fname), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["source", "target"])
            for pair in sorted(pairs):
                w.writerow(pair)
        lines.append("    " + yaml_quote(fname) + ":")
        lines.append("        type: [%s, %s, %s]" % (yaml_quote(et1), yaml_quote(rel), yaml_quote(et3)))
        lines.append("        cols: [0, 1]")
        lines.append("        skip: 1  # header row")
    lines += ["", "dynamics: {}  # TODO: add establishment and termination rules for the link types above (see README.md)", ""]
    with open(os.path.join(args.out, "model.yaml"), "w") as f:
        f.write("\n".join(lines))
    print("wrote %s/model.yaml: %d entities of %d types, %d relationship types, %d links of %d types"
          % (args.out, len(entities), len(et_entities), len(rel_label), sum(len(v) for v in links.values()), len(links)))
    if n_literal:
        print("  (%d triples with literal objects were ignored)" % n_literal)
    if n_self:
        print("  (%d self-links were ignored)" % n_self)


# ---------------------------------------------------------------------------------------------------------------
# tricl2rdf

def iri_local(label, used):
    """A safe, unique local name for a label."""
    base = re.sub(r"[^A-Za-z0-9_\-]+", "_", label).strip("_") or "x"
    if base[0].isdigit():
        base = "_" + base
    name = base
    k = 2
    while name in used and used[name] != label:
        name = "%s_%d" % (base, k)
        k += 1
    used[name] = label
    return name


def ttl_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def tricl2rdf(args):
    opener = gzip.open if args.input.endswith(".gz") else open
    with opener(args.input, "rb") as f:
        tree = ET.parse(f)
    root = tree.getroot()
    ns = {"g": root.tag.split("}")[0].strip("{")} if root.tag.startswith("{") else {"g": ""}

    def q(tag):
        return "{%s}%s" % (ns["g"], tag) if ns["g"] else tag

    used = {}
    out = ["@prefix : <%s> ." % args.base,
           "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .",
           "@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .",
           "@prefix tricl: <%s> ." % TRICL_NS, ""]
    node_iri = {}
    for node in root.iter(q("node")):
        label = node.get("label", node.get("id"))
        et = None
        for av in node.iter(q("attvalue")):
            if av.get("for") == "T":
                et = av.get("value")
        iri = ":" + iri_local(label, used)
        node_iri[node.get("id")] = iri
        line = iri
        if et:
            line += " a :" + iri_local(et, used)
        line += " ; rdfs:label " + ttl_string(label) + " ."
        out.append(line)
    out.append("")
    n_edges = 0
    for edge in root.iter(q("edge")):
        rel = None
        for av in edge.iter(q("attvalue")):
            if av.get("for") == "R":
                rel = av.get("value")
        if rel is None:
            continue
        s = node_iri.get(edge.get("source"))
        t = node_iri.get(edge.get("target"))
        if s is None or t is None:
            continue
        out.append("<< %s :%s %s >> tricl:from %s ; tricl:until %s ."
                   % (s, iri_local(rel, used), t, edge.get("start", "0"), edge.get("end", "")))
        n_edges += 1
    with open(args.out, "w") as f:
        f.write("\n".join(out) + "\n")
    print("wrote %s: %d nodes, %d link intervals" % (args.out, len(node_iri), n_edges))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    a = sub.add_parser("rdf2tricl", help="convert RDF data into a tricl config skeleton plus link csv files")
    a.add_argument("inputs", nargs="+")
    a.add_argument("--out", required=True)
    a.add_argument("--query", default=None)
    a.add_argument("--format", default=None)
    a.add_argument("--default-type", default="thing")
    a.add_argument("--type-map", action="append", default=[])
    a.add_argument("--relationship-map", action="append", default=[])
    a.add_argument("--symmetric", action="append", default=[])
    a.add_argument("--inverse", action="append", default=[])
    a.add_argument("--language", default="en")
    a.add_argument("--name", default="")
    a.set_defaults(func=rdf2tricl)
    b = sub.add_parser("tricl2rdf", help="convert a tricl gexf output file into Turtle-star with time intervals")
    b.add_argument("input")
    b.add_argument("--out", required=True)
    b.add_argument("--base", default="http://example.org/tricl/")
    b.set_defaults(func=tricl2rdf)
    args = ap.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
