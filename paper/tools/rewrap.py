"""rewrap.py -- one sentence per source line for LaTeX prose (the paper's source convention).

Adapted from the tool of the same name in the responsibility_theory project.

Usage: python3 paper/tools/rewrap.py FILE...           (rewrites the files in place)
       python3 paper/tools/rewrap.py --check FILE...   (lists files that violate the convention, exit code 1)

Figure and table captions are rewrapped too (one sentence per line inside \\caption{...}).

Only prose paragraphs are touched: consecutive non-blank lines outside environments whose
content is layout-sensitive (displays, tables, figures, listings, TikZ).  A paragraph is joined
and split again after each sentence end ('.', '?', '!' followed by white space and a capital or a
backslash command), never inside braces or inline math, and not after common abbreviations.
Lines with a comment, lines starting a new unit (\\item, \\paragraph, \\section, \\begin, ...) and
display math delimiters stay on their own lines.  Joining replaces newlines by single spaces,
which TeX reads the same, so the output PDF is unchanged (checked by the caller).
"""
import re
import sys

SKIP_ENVS = ("tabular", "tabular*", "longtable", "lstlisting", "verbatim", "tikzpicture",
             "align", "align*", "equation", "equation*", "gather", "gather*", "multline",
             "multline*", "figure", "figure*", "table", "table*", "CCSXML", "sideways",
             "minipage", "adjustbox")
ABBREV = ("e.g.", "i.e.", "cf.", "vs.", "Fig.", "Figs.", "Sect.", "Sects.", "App.", "Apps.",
          "Def.", "Defs.", "Eq.", "Eqs.", "al.", "resp.", "Thm.", "Prop.", "Cor.", "No.",
          "Lem.", "Ex.", "\\dots.", "approx.", "ca.", "w.r.t.", "p.", "pp.", "Vol.", "Ch.",
          "Sec.", "Secs.", "Ref.", "Refs.", "Eqs.", "Tab.", "init.", "est.")
# a line that must stay on its own: starts a unit or is a delimiter
OWN_LINE = re.compile(r"^\s*(\\item\b|\\(sub)*section\*?\b|\\paragraph\b|\\begin\b|\\end\b|"
                      r"\\\[|\\\]|\\caption\b|\\label\b|\\axhead\b|\\par\b|\\noindent\b|"
                      r"\\InputIfFileExists\b|\\input\b|\\clearpage\b|\\bigskip|\\medskip|"
                      r"\\smallskip|\\vspace|\\centering\b|\\mslabel|\\bibliograph|\\leanlist\b)")


# commands a prose line may start with; a line starting with any other command, a brace or a
# bracket (definitions, key lists, TikZ, ...) is not prose and stays as it is
INLINE = ("emph", "ms", "msphi", "ax", "axn", "lean", "leanloc", "leanfile", "lname", "op",
          "textit", "textbf", "textsc", "textsf", "texttt", "cite", "citet", "citep", "ref",
          "eqref", "nolinkurl", "url", "href", "Rb", "Rf", "srcHH", "srcNew", "srcMod", "srcHIH",
          "rw", "ding", "mbox", "enquote", "dots", "ldots", "S", "pageref", "hyperref", "Q", "Rl",
          "tricl", "ttt", "expit")
# macros for generated numbers (paper/generated/numbers.tex) are prose too
NUMBER_MACRO = re.compile(r"(Atop|Recovery|Macro|TestSigmoid)[A-Za-z]*$")


def is_prose(line):
    t = line.strip()
    if not t or t[0] in "{}[]&":
        return False
    m = re.match(r"\\([A-Za-z]+)", t)
    return m is None or m.group(1) in INLINE or NUMBER_MACRO.match(m.group(1)) is not None


def has_comment(line):
    return re.search(r"(?<!\\)%", line) is not None


def split_sentences(text):
    out, start, depth, math = [], 0, 0, False
    i = 0
    while i < len(text):
        c = text[i]
        if c == "\\":
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
        elif c == "$":
            math = not math
        elif c in ".?!" and depth == 0 and not math:
            m = re.match(r"[.?!]((?:''|\)|')*)( +)(?=[A-Z\\(`])", text[i:])
            if m:
                head = text[start:i + 1]
                if not any(head.endswith(a) for a in ABBREV) and \
                        not re.search(r"(^|[\s~(])[A-Z]\.$", head):   # initials
                    end = i + 1 + len(m.group(1))
                    out.append(text[start:end])
                    start = end + len(m.group(2))
                    i = start
                    continue
        i += 1
    out.append(text[start:])
    return [s for s in out if s.strip()]


def rewrap(lines):
    out, para, env_stack = [], [], []

    def flush():
        if not para:
            return
        indent = re.match(r"\s*", para[0]).group(0)
        cont = re.match(r"\s*", para[1]).group(0) if len(para) > 1 else indent
        text = " ".join(l.strip() for l in para)
        sents = split_sentences(text)
        out.append(indent + sents[0])
        out.extend(cont + s for s in sents[1:])
        para.clear()

    hold = 0                       # open braces of a kept non-prose line (definitions, key lists)

    def net(l):
        l = re.sub(r"\\[{}]|(?<!\\)%.*", "", l)
        return l.count("{") - l.count("}")

    for line in lines:
        if hold > 0:
            flush()
            out.append(line)
            hold += net(line)
            continue
        begins = re.findall(r"\\begin\{([^}]*)\}", line)
        ends = re.findall(r"\\end\{([^}]*)\}", line)
        in_skip = any(e in SKIP_ENVS for e in env_stack) or "\\[" in "".join(env_stack)
        opens_skip = any(b in SKIP_ENVS for b in begins) or (
            re.search(r"\\\[", line) and not re.search(r"\\\]", line))
        if in_skip or opens_skip or not line.strip() or has_comment(line) or OWN_LINE.match(line) \
                or re.search(r"\\\[|\\\]|\\\\\s*$", line) or (not para and not is_prose(line)) \
                or (para and not is_prose(line) and not re.match(r"\s", line)):
            flush()
            out.append(line)
            if not in_skip and not opens_skip:
                hold = max(0, net(line))
        else:
            para.append(line)
        for b in begins:
            env_stack.append(b)
        for e in ends:
            if e in env_stack:
                env_stack.reverse(); env_stack.remove(e); env_stack.reverse()
        if re.search(r"\\\[", line) and not re.search(r"\\\]", line):
            env_stack.append("\\[")
        elif re.search(r"\\\]", line) and "\\[" in env_stack and not re.search(r"\\\[", line):
            env_stack.remove("\\[")
        # a unit-starting line such as \item text... continues as a paragraph
        if OWN_LINE.match(line) and re.match(r"^\s*\\(item|paragraph|axhead)\b", line) \
                and not has_comment(line) and not opens_skip and not in_skip:
            out.pop()
            para.append(line)
    flush()
    return out


def rewrap_captions(src):
    """Put each sentence of every \\caption{...} on its own line."""
    out, pos = [], 0
    for m in re.finditer(r"\\caption\{", src):
        if m.start() < pos:
            continue
        i, depth = m.end(), 1
        while depth:
            c = src[i]
            if c == "\\":
                i += 2
                continue
            depth += (c == "{") - (c == "}")
            i += 1
        body = src[m.end():i - 1]
        if re.search(r"(?<!\\)%", body):        # leave captions with comments alone
            continue
        text = " ".join(l.strip() for l in body.split("\n"))
        out.append(src[pos:m.end()] + "\n".join(split_sentences(text)) + "}")
        pos = i
    out.append(src[pos:])
    return "".join(out)


def process(src):
    trailing = src.endswith("\n")
    new = "\n".join(rewrap(src.split("\n")[:-1] if trailing else src.split("\n")))
    return rewrap_captions(new + ("\n" if trailing else ""))


if __name__ == "__main__":
    check = "--check" in sys.argv[1:]
    bad = []
    for path in [a for a in sys.argv[1:] if a != "--check"]:
        src = open(path, encoding="utf-8").read()
        new = process(src)
        if check:
            if new != src:
                bad.append(path)
        elif new != src:
            open(path, "w", encoding="utf-8").write(new)
    if check:
        for path in bad:
            print("not one sentence per line:", path)
        sys.exit(1 if bad else 0)
