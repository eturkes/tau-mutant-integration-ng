#!/usr/bin/env python3
"""build_synthesis_model_figure.py -- schematic for the synthesis.Rmd lead model.

Render the "Architect-and-Modifier" unified holistic model (M-UNIFIED) as a
standalone wiring diagram: amyloid as the ARCHITECT that templates the
homeostatic->DAM activation trajectory, tau as a MODIFIER firing at a
supra-additive GSK3b coincidence (AND-gate) node that drives two
post-transcriptional channels -- a within-state progression RATE synergy and a
WIRING cascade (GSK3b -| Mapk14 -> Myc/Foxo3; the amyloid NF-kB increment
attenuated) -- neither of which rewrites the DAM differentially-expressed gene
set, whose intercellular output stays TREM2/APP-fragment clearance.

PURE presentation: the node set, layout and signed edges are the authored model
topology (curated + critique-corrected output of the holistic-models workflow:
all five interaction TFs originate from a SINGLE GSK3b seed in CARNIVAL; the
NF-kB edge and the p38->GSK3b feedback are recovered OFF-interaction and are
flagged PROPOSED). Edge styling encodes sign (activation arrowhead / inhibition
bar), proposed status (dashed), and temporal order (dotted). Mirrors
build_capstone_figure.py (matplotlib Agg, DejaVu/Arial, 300 dpi, tracked PNG).

Output: storage/figures/synthesis_model.png  (tracked; re-run to refresh).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Patch
from matplotlib.lines import Line2D

OUT = "storage/figures/synthesis_model.png"

# ---- authored model topology ----------------------------------------------
# id -> (x, y, label, kind)
NODES = {
    "amyloid":     (0.8, 8.4, "Amyloid\n(App NL-G-F)",            "driver"),
    "tau":         (0.8, 5.0, "Tau\n(P301S)",                     "driver"),
    "gsk3b":       (3.6, 5.0, "GSK3β\ncoincidence AND-gate",  "gate"),
    "dam_traj":    (3.6, 8.4, "Homeostatic→DAM\nactivation trajectory", "state"),
    "dam_prog":    (6.5, 8.4, "DAM gene\nprogramme",              "state"),
    "clearance":   (9.4, 8.4, "TREM2 / APP-fragment\nclearance output", "output"),
    "progression": (6.5, 5.0, "Progression\nRATE ↑",          "state"),
    "nfkb":        (9.4, 6.5, "NF-κB increment\nATTENUATED",   "signal"),
    "mapk14":      (3.6, 2.2, "p38 / Mapk14",                     "signal"),
    "myc":         (6.5, 2.2, "Myc / Foxo3\nbiosynthesis ↓",  "signal"),
}

# (src, dst, kind, proposed, label, rad)   kind in {act, inh, prec, fb}
EDGES = [
    ("amyloid", "dam_traj",    "act", False, "architect  +7.9 / +10.3", 0.0),
    ("dam_traj", "dam_prog",   "act", False, "Cst7/Itgax/Lpl ↑",   0.0),
    ("dam_prog", "clearance",  "act", False, "Apoe-Trem2 · App-Cd74", 0.0),
    ("amyloid", "gsk3b",       "act", False, "input 1 (alone −)",  0.0),
    ("tau", "gsk3b",           "act", False, "input 2 (alone −)",  0.0),
    ("gsk3b", "progression",   "act", False, "channel 1 · rate +2.46", 0.0),
    ("progression", "dam_traj","act", False, "faster traversal",        0.30),
    ("gsk3b", "mapk14",        "inh", False, "channel 2 · wiring",  0.0),
    ("mapk14", "myc",          "act", False, "p38 normally activates → net ↓", 0.0),
    ("dam_prog", "nfkb",       "act", False, "amyloid increment",       0.0),
    ("gsk3b", "nfkb",          "inh", True,  "attenuated (proposed)",  -0.34),
    ("amyloid", "tau",         "prec", True, "inferred order",          0.0),
    ("mapk14", "gsk3b",        "fb",  True,  "proposed feedback",       0.40),
]

COL = {"driver": "#f6c89a", "gate": "#ffd54a", "state": "#bcd9ee",
       "signal": "#d9c2e9", "output": "#bfe3c6"}
KIND_LABEL = {"driver": "perturbation", "gate": "coincidence node",
              "state": "cell-state axis", "signal": "signalling node",
              "output": "intercellular output"}
ECOL = {"act": "#333333", "inh": "#b2182b", "prec": "#7a7a7a", "fb": "#6a51a3"}

NW, NH = 2.04, 1.0  # node box width / height

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.8,
})
fig, ax = plt.subplots(figsize=(13.6, 8.8))
ax.set_xlim(-0.35, 10.95)
ax.set_ylim(1.15, 9.45)
ax.axis("off")

# ---- nodes -----------------------------------------------------------------
patches = {}
for nid, (x, y, label, kind) in NODES.items():
    p = FancyBboxPatch((x - NW / 2, y - NH / 2), NW, NH,
                       boxstyle="round,pad=0.02,rounding_size=0.14",
                       linewidth=1.4, edgecolor="#2a2a2a",
                       facecolor=COL[kind], zorder=3)
    ax.add_patch(p)
    patches[nid] = p
    fw = "bold" if kind == "gate" else "normal"
    ax.text(x, y, label, ha="center", va="center", fontsize=9.0,
            fontweight=fw, color="#161616", zorder=4, linespacing=1.05)

# ---- edges -----------------------------------------------------------------
def edge_style(kind, proposed):
    if kind == "inh":
        astyle = "|-|,widthA=0.0,widthB=0.55"
    else:
        astyle = "-|>"
    ls = ":" if kind == "prec" else ("--" if proposed else "-")
    alpha = 0.80 if proposed else 0.95
    lw = 1.7 if proposed else 2.1
    return astyle, ls, alpha, lw

for src, dst, kind, proposed, label, rad in EDGES:
    xa, ya = NODES[src][0], NODES[src][1]
    xb, yb = NODES[dst][0], NODES[dst][1]
    astyle, ls, alpha, lw = edge_style(kind, proposed)
    arr = FancyArrowPatch((xa, ya), (xb, yb), patchA=patches[src],
                          patchB=patches[dst], shrinkA=3, shrinkB=3,
                          connectionstyle=f"arc3,rad={rad}", arrowstyle=astyle,
                          mutation_scale=20, linewidth=lw, linestyle=ls,
                          color=ECOL[kind], alpha=alpha, zorder=2,
                          joinstyle="round", capstyle="round")
    ax.add_patch(arr)
    # label at the (possibly bowed) midpoint
    mx, my = (xa + xb) / 2.0, (ya + yb) / 2.0
    dx, dy = xb - xa, yb - ya
    L = max(np.hypot(dx, dy), 1e-6)
    nx, ny = -dy / L, dx / L              # unit perpendicular
    off = 0.30 + 0.5 * rad * L            # base offset + arc bulge
    lx, ly = mx + nx * off, my + ny * off
    ax.text(lx, ly, label, ha="center", va="center", fontsize=7.1,
            color=ECOL[kind], zorder=5, fontstyle="italic",
            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none", alpha=0.78))

# annotation: the single-seed TF fan-out (corrected CARNIVAL topology)
ax.text(6.5, 1.32,
        "all five interaction TFs originate from one GSK3β seed (CARNIVAL): "
        "Myc/Foxo3 via p38; Clock/Sfpq direct; Tbp via Csnk2b. "
        "Coordinated suppression also of Creb1/Tp53/Jun.",
        ha="center", va="center", fontsize=7.2, color="#444444", fontstyle="italic")

# ---- title + legend --------------------------------------------------------
fig.suptitle("The lead holistic model: amyloid as architect, tau as a rate-and-wiring modifier "
             "at a supra-additive GSK3β node",
             x=0.5, y=0.975, fontsize=12.6, fontweight="bold")
fig.text(0.5, 0.937,
         "Tau acts only on the amyloid-built trajectory: a GSK3β AND-gate (suppressed by each insult alone, fired by the combination) drives a "
         "progression-RATE channel and a post-translational WIRING channel, neither rewrites the DAM gene set, so the interaction is "
         "static-DE-null yet kinase/causal/trajectory-positive.",
         ha="center", fontsize=8.5, color="#444444")

node_handles = [Patch(facecolor=COL[k], edgecolor="#2a2a2a", label=KIND_LABEL[k])
                for k in ("driver", "gate", "state", "signal", "output")]
edge_handles = [
    Line2D([0], [0], color=ECOL["act"], lw=2.1, marker=">", markersize=7,
           markerfacecolor=ECOL["act"], label="activation / promotes"),
    Line2D([0], [0], color=ECOL["inh"], lw=2.1, marker="|", markersize=10,
           markeredgewidth=2.2, label="inhibition / suppression"),
    Line2D([0], [0], color=ECOL["prec"], lw=1.8, ls=":", label="inferred temporal order"),
    Line2D([0], [0], color="#555555", lw=1.7, ls="--", label="proposed (off-interaction / reach)"),
]
leg1 = ax.legend(handles=node_handles, loc="lower left", bbox_to_anchor=(0.0, 0.0),
                 frameon=False, fontsize=7.8, handlelength=1.3, title="Node type",
                 title_fontsize=8.2, labelspacing=0.5)
leg1.get_title().set_fontweight("bold")
ax.add_artist(leg1)
leg2 = ax.legend(handles=edge_handles, loc="lower right", bbox_to_anchor=(1.0, 0.0),
                 frameon=False, fontsize=7.8, handlelength=2.0, title="Edge type",
                 title_fontsize=8.2, labelspacing=0.5)
leg2.get_title().set_fontweight("bold")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=300, facecolor="white", bbox_inches="tight")
print(f"[synthesis-model-fig] wrote {OUT}  ({fig.get_size_inches()[0]:.1f}x{fig.get_size_inches()[1]:.1f} in @300dpi)")
