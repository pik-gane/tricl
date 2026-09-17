#!/usr/bin/env python3
"""Render a movie of a temporal network written by tricl (``--entities-out`` and ``--links-out``).

Usage:
    python3 python/tricl_movie.py ENTITIES.csv LINKS.csv --out MOVIE.mp4 [options]

The links file is the interval list written by ``tricl --links-out``: one row per directed link with columns
``source,relationship,target,start,end``. The entities file (``tricl --entities-out``) has columns ``id,label,type``.

Frames are taken at equally spaced model times between ``--t0`` and ``--t1`` (default: the whole run). At each frame,
the links with ``start <= t < end`` are drawn (links still existing at the end of the run are drawn until the end).
Two links in opposite directions (as tricl writes for symmetric relationship types) are drawn as one line, a single
directed link as an arrow.

Node positions come from a Fruchterman-Reingold layout of the *time-aggregated* graph, in which every pair of
entities is attracted in proportion to the total time it was linked. This keeps the picture still, so that the
movement one sees is the dynamics of the links and states, not the layout algorithm. ``--layout dynamic`` instead
recomputes the layout for every frame, warm-started from the previous one, for slowly changing networks.

Nodes are coloured by entity type, or with ``--state R [R ...]`` by their current link of relationship type R:
tricl models often encode an entity's state as a link to a "hub" entity (e.g. ``[agent, is, active]`` or
``[agent, is infectious for, covid-19]``), and this option colours each source entity by that link (the state is
named after the target entity, the relationship type, or both, whichever distinguishes the states). The hub
entities themselves are then not drawn.

The output is an mp4/webm/mkv movie (needs the ``ffmpeg`` program, or the ``imageio-ffmpeg`` package which bundles
it), an animated gif (via Pillow), or a directory of png frames. Needs numpy and matplotlib.

Examples:
    tricl config_files/granovetter_helfmann.yaml --seed 1 --E 5000 --entities-out entities.csv --links-out links.csv
    python3 python/tricl_movie.py entities.csv links.csv --state is "is not" --out granovetter.mp4

    tricl config_files/sir_sd.yaml --seed 1 --entities-out entities.csv --links-out links.csv
    python3 python/tricl_movie.py entities.csv links.csv --draw "regularly meets" \\
        --state "is susceptible to" "is infectious for" "has recovered from" "has died of" --out sir.gif

The interval list itself is the common input format of temporal network libraries (each row is one edge with a
start and end time, in a plain csv table); see the README.
"""

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

# The colours follow a validated categorical palette (fixed order, colour-vision-deficiency safe for adjacent slots);
# the "dark" steps are the same hues adjusted for a dark surface.
PALETTES = {
    "light": {
        "series": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"],
        "surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781", "grid": "#e1e0d9",
    },
    "dark": {
        "series": ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"],
        "surface": "#1a1a19", "ink": "#ffffff", "ink2": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a",
    },
}


class Network:
    """The entities and the interval list, as arrays."""

    def __init__(self, entities_file, links_file):
        self.labels, self.types = [], []
        with open(entities_file, newline="") as f:
            for row in csv.DictReader(f):
                self.labels.append(row["label"])
                self.types.append(row["type"])
        self.index = {label: i for i, label in enumerate(self.labels)}
        if len(self.index) != len(self.labels):
            raise ValueError("entity labels are not unique")
        src, dst, rel, start, end = [], [], [], [], []
        self.relationships = []
        rel_index = {}
        with open(links_file, newline="") as f:
            for row in csv.DictReader(f):
                try:
                    s, d = self.index[row["source"]], self.index[row["target"]]
                except KeyError as e:
                    raise ValueError("unknown entity %s in %s" % (e, links_file))
                r = rel_index.get(row["relationship"])
                if r is None:
                    r = rel_index[row["relationship"]] = len(self.relationships)
                    self.relationships.append(row["relationship"])
                src.append(s); dst.append(d); rel.append(r)
                start.append(float(row["start"])); end.append(float(row["end"]))
        self.src, self.dst, self.rel = (np.array(a, dtype=np.int64) for a in (src, dst, rel))
        self.start, self.end = np.array(start, dtype=float), np.array(end, dtype=float)
        self.t_end = float(self.end.max()) if len(self.end) else 0.0
        self.rel_index = rel_index

    @property
    def n(self):
        return len(self.labels)

    def alive(self, t):
        """Boolean mask of the links existing at time t (links ending at the end of the run count as still existing)."""
        return (self.start <= t) & ((self.end > t) | (self.end >= self.t_end))


def fr_layout(n, u, v, w, pos=None, iterations=50, temperature=None, seed=0):
    """Fruchterman-Reingold layout (numpy, O(n^2) per iteration).

    u, v, w: edge endpoints and weights (attraction is proportional to the weight). pos: starting positions
    (random if None). The positions are not normalised (see normalize()).
    """
    rng = np.random.default_rng(seed)
    if pos is None:
        pos = rng.uniform(-1, 1, size=(n, 2))
    else:
        pos = np.array(pos, dtype=float)
    if n <= 1:
        return np.zeros((max(n, 0), 2))
    k = 2.0 / math.sqrt(n)  # ideal edge length in a box of side 2
    t = 0.2 if temperature is None else temperature
    dt = t / (iterations + 1)
    w = np.asarray(w, dtype=float)
    for _ in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]
        dist2 = (delta ** 2).sum(axis=-1)
        np.fill_diagonal(dist2, 1.0)
        dist2 = np.maximum(dist2, 1e-4)
        disp = (delta * (k * k / dist2)[:, :, None]).sum(axis=1)  # repulsion from all other nodes
        if len(u):
            d = pos[u] - pos[v]
            dl = np.maximum(np.sqrt((d ** 2).sum(axis=1)), 1e-6)
            f = (dl * w / k)[:, None] * d  # attraction along edges (force dl^2/k times weight)
            np.add.at(disp, u, -f)
            np.add.at(disp, v, f)
        length = np.maximum(np.sqrt((disp ** 2).sum(axis=1)), 1e-9)
        pos = pos + disp / length[:, None] * np.minimum(length, t)[:, None]
        t -= dt
    return pos


def normalize(pos):
    """Centre the positions and scale them so that 98% of them lie within the box [-1, 1]^2 (the rest is clipped)."""
    pos = pos - np.nanmean(pos, axis=0)
    scale = np.nanpercentile(np.abs(pos), 98)
    if scale > 0:
        pos = pos / scale
    return np.clip(pos, -1.08, 1.08)


def aggregate_edges(net, mask, t0, t1, weighted=True):
    """Unordered pairs linked during [t0, t1] among the rows selected by mask, with total link duration as weight."""
    rows = np.flatnonzero(mask)
    if not len(rows):
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0)
    a = np.minimum(net.src[rows], net.dst[rows])
    b = np.maximum(net.src[rows], net.dst[rows])
    if weighted:
        overlap = np.maximum(0.0, np.minimum(net.end[rows], t1) - np.maximum(net.start[rows], t0))
    else:
        overlap = np.ones(len(rows))
    key = a * net.n + b
    keys, inverse = np.unique(key, return_inverse=True)
    weight = np.bincount(inverse, weights=overlap, minlength=len(keys))
    keep = weight > 0
    keys, weight = keys[keep], weight[keep]
    return keys // net.n, keys % net.n, weight / weight.max()


def find_ffmpeg():
    exe = shutil.which("ffmpeg")
    if exe:
        return exe
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split("Examples:")[1] if "Examples:" in __doc__ else None)
    ap.add_argument("entities", help="csv file written by tricl --entities-out")
    ap.add_argument("links", help="csv file written by tricl --links-out")
    ap.add_argument("--out", default="frames", help="output: FILE.mp4/.webm/.mkv (ffmpeg), FILE.gif (Pillow), or a directory for png frames (default: frames)")
    ap.add_argument("--frames", type=int, default=200, help="number of frames (default 200)")
    ap.add_argument("--dt", type=float, default=None, help="model time between frames (alternative to --frames)")
    ap.add_argument("--fps", type=float, default=20, help="frames per second (default 20)")
    ap.add_argument("--t0", type=float, default=0.0, help="first frame time (default 0)")
    ap.add_argument("--t1", type=float, default=None, help="last frame time (default: end of the run)")
    ap.add_argument("--layout", choices=["static", "dynamic"], default="static", help="static: one layout of the time-aggregated graph (default); dynamic: relayout every frame")
    ap.add_argument("--iterations", type=int, default=None, help="layout iterations (default 100 static, 5 per frame dynamic)")
    ap.add_argument("--layout-relationships", nargs="+", default=None, metavar="R", help="relationship types that hold the layout together (default: the drawn ones)")
    ap.add_argument("--draw", nargs="+", default=None, metavar="R", help="relationship types to draw (default: all except state relationships)")
    ap.add_argument("--state", nargs="+", default=None, metavar="R", help="colour entities by their current link of these relationship types (see above)")
    ap.add_argument("--hide", nargs="+", default=[], metavar="TYPE", help="entity types not to draw")
    ap.add_argument("--labels", action="store_true", help="write the entity labels next to the nodes")
    ap.add_argument("--title", default=None, help="title written above the picture")
    ap.add_argument("--size", type=float, nargs=2, default=(8, 6), metavar=("W", "H"), help="figure size in inches (default 8 6)")
    ap.add_argument("--dpi", type=int, default=100, help="dots per inch (default 100)")
    ap.add_argument("--node-size", type=float, default=None, help="node diameter in pixels (default: between 8 and 24, depending on the number of nodes)")
    ap.add_argument("--edge-width", type=float, default=None, help="edge width in pixels (default: 2 for few edges, thinner for many)")
    ap.add_argument("--dark", action="store_true", help="dark surface")
    ap.add_argument("--seed", type=int, default=0, help="random seed of the layout")
    ap.add_argument("--keep-frames", default=None, metavar="DIR", help="also keep the png frames in this directory")
    args = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    pal = PALETTES["dark" if args.dark else "light"]
    net = Network(args.entities, args.links)
    t0 = args.t0
    t1 = net.t_end if args.t1 is None else args.t1
    if not (t1 >= t0):
        raise SystemExit("empty time window")
    if args.dt:
        times = np.arange(t0, t1 + 1e-9 * max(1.0, abs(t1)), args.dt)
    else:
        times = np.linspace(t0, t1, max(args.frames, 1))

    def rel_ids(names, what):
        ids = []
        for name in names:
            if name not in net.rel_index:  # (a relationship type may exist in the model but have no links in this run)
                print("%s: no links of relationship type '%s' in %s (types found: %s)" % (what, name, args.links, ", ".join(net.relationships)))
                continue
            ids.append(net.rel_index[name])
        return ids

    # states:
    state_rels = rel_ids(args.state or [], "--state")
    hidden = set(args.hide)
    if state_rels:
        state_mask = np.isin(net.rel, state_rels)
        for e in np.unique(net.dst[state_mask]):
            hidden.add(net.types[e])  # hub entities are not drawn
        combos = sorted({(int(r), int(d)) for r, d in zip(net.rel[state_mask], net.dst[state_mask])}, key=lambda c: (state_rels.index(c[0]), c[1]))
        many_rels = len({r for r, d in combos}) > 1
        many_targets = len({d for r, d in combos}) > 1
        combo_index = {c: i for i, c in enumerate(combos)}
        state_names = [((net.relationships[r] + " ") if many_rels or not many_targets else "") + (net.labels[d] if many_targets else "") for r, d in combos]
        state_names = [s.strip() for s in state_names]
        combo_of_row = np.full(len(net.rel), -1, dtype=np.int64)
        for i in np.flatnonzero(state_mask):
            combo_of_row[i] = combo_index[(int(net.rel[i]), int(net.dst[i]))]
    else:
        state_names, combo_of_row = [], None

    shown = np.array([t not in hidden for t in net.types])
    if not shown.any():
        raise SystemExit("no entities left to draw")
    drawn_types = []
    for t in net.types:
        if t not in hidden and t not in drawn_types:
            drawn_types.append(t)

    # relationship types to draw:
    if args.draw is not None:
        draw_rels = rel_ids(args.draw, "--draw")
    else:
        draw_rels = [r for r in range(len(net.relationships)) if r not in state_rels]
    drawable = np.isin(net.rel, draw_rels) & shown[net.src] & shown[net.dst]
    draw_rels = [r for r in draw_rels if (net.rel[drawable] == r).any()]
    layout_rels = rel_ids(args.layout_relationships, "--layout-relationships") if args.layout_relationships else draw_rels
    layout_mask = np.isin(net.rel, layout_rels) & shown[net.src] & shown[net.dst]

    # colours: node categories first, then relationship types, in the fixed palette order:
    node_cats = state_names if state_rels else drawn_types
    slots = list(pal["series"])
    node_colors = {}
    for c in node_cats:
        node_colors[c] = slots.pop(0) if slots else pal["muted"]
    edge_colors = {}
    if len(draw_rels) > 1:
        for r in draw_rels:
            edge_colors[r] = slots.pop(0) if slots else pal["muted"]
    else:
        for r in draw_rels:
            edge_colors[r] = pal["muted"]

    # layout:
    iterations = args.iterations or (100 if args.layout == "static" else 5)
    u, v, w = aggregate_edges(net, layout_mask, t0, t1)
    print("layout of %d shown entities and %d aggregated pairs ..." % (shown.sum(), len(u)), flush=True)
    pos = fr_layout(net.n, u, v, w, iterations=iterations if args.layout == "static" else 100, seed=args.seed)
    pos[~shown] = np.nan
    pos = normalize(pos)

    # figure:
    W, H = args.size
    fig = plt.figure(figsize=(W, H), dpi=args.dpi, facecolor=pal["surface"])
    legend_w = 0.24 if (len(node_cats) > 1 or len(draw_rels) > 1) else 0.02
    ax = fig.add_axes([0.02, 0.06, 0.96 - legend_w, 0.86])
    ax.set_facecolor(pal["surface"])
    ax.set_xlim(-1.12, 1.12); ax.set_ylim(-1.12, 1.12)
    ax.set_aspect("equal"); ax.axis("off")
    px = 72.0 / args.dpi  # points per pixel
    n_shown = int(shown.sum())
    node_px = args.node_size or min(24.0, max(8.0, 300.0 / math.sqrt(n_shown)))
    n_edges_typical = max(1, int(np.median([net.alive(t)[drawable].sum() for t in times[:: max(1, len(times) // 10)]])))
    edge_px = args.edge_width or (2.0 if n_edges_typical <= 300 else max(0.6, 2.0 * math.sqrt(300.0 / n_edges_typical)))
    edge_alpha = 0.7 if n_edges_typical <= 300 else 0.45
    lines = {r: LineCollection([], linewidths=edge_px * px, colors=edge_colors[r], alpha=edge_alpha, capstyle="round", zorder=1) for r in draw_rels}
    for lc in lines.values():
        ax.add_collection(lc)
    arrows = []
    nodes = ax.scatter([], [], s=(node_px * px) ** 2, edgecolors=pal["surface"], linewidths=2 * px, zorder=3)
    label_artists = []
    if args.labels:
        for i in np.flatnonzero(shown):
            label_artists.append(ax.annotate(net.labels[i], (0, 0), xytext=(node_px / 2 + 2, 0), textcoords="offset pixels",
                                             fontsize=7, color=pal["ink2"], va="center", zorder=4))
    if args.title:
        fig.text(0.02, 0.965, args.title, fontsize=13, color=pal["ink"], va="center", fontfamily="sans-serif")
    subtitle = fig.text(0.02, 0.925 if args.title else 0.955, "", fontsize=10, color=pal["ink2"], va="center")
    # timeline: a hairline with the elapsed part in secondary ink
    fig.add_artist(Line2D([0.02, 0.98 - legend_w], [0.03, 0.03], color=pal["grid"], linewidth=1, transform=fig.transFigure))
    progress = Line2D([0.02, 0.02], [0.03, 0.03], color=pal["ink2"], linewidth=2, transform=fig.transFigure, solid_capstyle="round")
    fig.add_artist(progress)
    # legend (always present for two or more categories; text in ink, never in the data colour):
    handles = []
    if len(node_cats) > 1:
        for c in node_cats:
            handles.append(Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=node_colors[c], markeredgecolor=pal["surface"], label=c))
    if len(draw_rels) > 1:
        for r in draw_rels:
            handles.append(Line2D([], [], color=edge_colors[r], linewidth=2, label=net.relationships[r]))
    if handles:
        leg = fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.98 - legend_w, 0.5), frameon=False, fontsize=9,
                         labelcolor=pal["ink"], handletextpad=0.6, borderaxespad=0)
        leg.set_in_layout(False)

    # frames:
    frames_dir = args.keep_frames or (args.out if not os.path.splitext(args.out)[1] else tempfile.mkdtemp(prefix="tricl_frames_"))
    os.makedirs(frames_dir, exist_ok=True)
    frame_files = []
    prev_pos = pos.copy()
    for i, t in enumerate(times):
        alive = net.alive(t)
        if args.layout == "dynamic" and i > 0:
            u, v, w = aggregate_edges(net, alive & layout_mask, t, t, weighted=False)
            filled = np.where(np.isnan(prev_pos), 0.0, prev_pos)
            new_pos = fr_layout(net.n, u, v, w, pos=filled, iterations=iterations, temperature=0.03, seed=args.seed)
            new_pos[~shown] = np.nan
            pos = normalize(new_pos)
            prev_pos = pos.copy()
        # edges: mutual links as lines, one-way links as arrows
        rows = np.flatnonzero(alive & drawable)
        a, b, r = net.src[rows], net.dst[rows], net.rel[rows]
        key = (r * net.n + a) * net.n + b
        rkey = (r * net.n + b) * net.n + a
        mutual = np.isin(rkey, key)
        for art in arrows:
            art.remove()
        arrows = []
        for rr in draw_rels:
            sel = (r == rr) & mutual & (a < b)
            segs = np.stack([pos[a[sel]], pos[b[sel]]], axis=1) if sel.any() else np.zeros((0, 2, 2))
            lines[rr].set_segments(segs)
            sel = (r == rr) & ~mutual
            if sel.any():
                p, q = pos[a[sel]], pos[b[sel]]
                d = q - p
                dl = np.maximum(np.sqrt((d ** 2).sum(axis=1)), 1e-9)[:, None]
                shrink = node_px / 2 * px / 72 / (W * 0.96 - legend_w * W) * 2.24  # data units per node radius, roughly
                q2 = q - d / dl * shrink
                arrows.append(ax.quiver(p[:, 0], p[:, 1], (q2 - p)[:, 0], (q2 - p)[:, 1], angles="xy", scale_units="xy", scale=1,
                                        width=edge_px * 0.0025, headwidth=5, headlength=6, headaxislength=5, minlength=0,
                                        color=edge_colors[rr], alpha=edge_alpha, zorder=2))
        # nodes:
        if state_rels:
            cat = np.full(net.n, -1, dtype=np.int64)
            srows = np.flatnonzero(alive & (combo_of_row >= 0))
            srows = srows[np.argsort(net.start[srows], kind="stable")]  # the latest state link wins
            cat[net.src[srows]] = combo_of_row[srows]
            colors = [node_colors[state_names[c]] if c >= 0 else pal["muted"] for c in cat[shown]]
        else:
            colors = [node_colors[net.types[j]] for j in np.flatnonzero(shown)]
        nodes.set_offsets(pos[shown])
        nodes.set_facecolors(colors)
        for art, j in zip(label_artists, np.flatnonzero(shown)):
            art.xy = pos[j]
        subtitle.set_text("t = %.6g    %d links" % (t, len(rows) - mutual.sum() // 2))
        progress.set_xdata([0.02, 0.02 + (0.96 - legend_w) * ((t - t0) / (t1 - t0) if t1 > t0 else 1.0)])
        fn = os.path.join(frames_dir, "frame_%05d.png" % i)
        fig.savefig(fn, dpi=args.dpi, facecolor=pal["surface"])
        frame_files.append(fn)
        if (i + 1) % 20 == 0 or i + 1 == len(times):
            print("  %d/%d frames" % (i + 1, len(times)), flush=True)
    plt.close(fig)

    # assemble:
    ext = os.path.splitext(args.out)[1].lower()
    if ext == "":
        print("wrote %d frames to %s" % (len(frame_files), frames_dir))
    elif ext == ".gif":
        from PIL import Image
        images = [Image.open(fn).convert("P", palette=Image.ADAPTIVE) for fn in frame_files]
        images[0].save(args.out, save_all=True, append_images=images[1:], duration=int(round(1000 / args.fps)), loop=0)
        print("wrote", args.out)
    else:
        ffmpeg = find_ffmpeg()
        if not ffmpeg:
            raise SystemExit("ffmpeg not found (install it, or 'pip install imageio-ffmpeg'); frames are in %s, or use --out FILE.gif" % frames_dir)
        codec = ["-c:v", "libvpx-vp9", "-b:v", "0", "-crf", "30"] if ext == ".webm" else ["-c:v", "libx264", "-crf", "20"]
        cmd = [ffmpeg, "-y", "-loglevel", "error", "-framerate", str(args.fps), "-i", os.path.join(frames_dir, "frame_%05d.png"),
               *codec, "-pix_fmt", "yuv420p", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", args.out]
        subprocess.run(cmd, check=True)
        print("wrote", args.out)
    if ext != "" and not args.keep_frames:
        shutil.rmtree(frames_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
