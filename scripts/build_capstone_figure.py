#!/usr/bin/env python3
"""build_capstone_figure.py -- arc P capstone deliverable C3.

Render the project-closing convergence figure as a polished standalone PNG from
the two C2 TSVs (storage/results/capstone_{convergence_matrix,contest_summary}.tsv).
PURE presentation: reads the curated, audit-sourced matrix; computes nothing.

Two panels:
  A  Landscape view -- a convergence heatmap, 12 mouse evidence-class rows x 3
     axis columns, with the single Human row held in its own separated band
     below a divider (evidence-class-incommensurable per H7; never counted in
     the contests). Cells are coloured by sign x grade (diverging - blue / 0
     grey / + red, intensity = grade) with the cell token printed IN each cell
     as a redundant, grayscale- and colourblind-safe encoding. Meta cells:
     purple = mixed (sign-divergent), grey = null (tested, no effect), gold+hatch
     = reframes, near-white = n.a. (axis not probed).
  B  Contest view -- the three ledger-adjudicated margins (sealed sec-17), each
     bar annotated with the winning hypothesis and its Strong/Moderate/Suggestive
     support counts.

Output: storage/figures/capstone_convergence.png (tracked; not regenerated
in-knit -- re-run this script to refresh).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import pandas as pd

RES = "storage/results"
OUT = "storage/figures/capstone_convergence.png"

# ---- load the curated synthesis -------------------------------------------
mat = pd.read_csv(os.path.join(RES, "capstone_convergence_matrix.tsv"), sep="\t")
con = pd.read_csv(os.path.join(RES, "capstone_contest_summary.tsv"), sep="\t")

AXES = ["amyloid_activation", "synaptic_suppression", "interaction_metabolic"]
AXIS_LABELS = ["Amyloid\nactivation", "Synaptic\nsuppression", "Tau×amyloid\ninteraction"]

# narrative-arc row order (D->O); base first, Human held out for its own band
MOUSE_ROWS = [
    ("Expression DE (4 modalities)", "base"),
    ("Pathway / module", "D"),
    ("TF activity", "E"),
    ("Kinase activity", "F"),
    ("CCC / ligand-receptor", "G"),
    ("NF-kB attenuation", "I"),
    ("Causal topology", "J"),
    ("SCENIC regulons", "K"),
    ("Composition", "L"),
    ("Dynamics / progression", "M"),
    ("Specificity", "N"),
    ("Dynamics-DE", "O"),
]
HUMAN_ROW = ("Human cross-species", "rmd 18")

# integrity: the figure must cover exactly the matrix's evidence classes
got = set(mat["evidence_class"])
want = {r[0] for r in MOUSE_ROWS} | {HUMAN_ROW[0]}
assert got == want, f"row-set mismatch: missing {want - got}, extra {got - want}"

cell_of = {(r.evidence_class, r.axis): r for r in mat.itertuples()}

# ---- sign x grade -> colour ------------------------------------------------
POS = {"Strong": "#b2182b", "Moderate": "#ef8a62", "Suggestive": "#fddbc7"}
NEG = {"Strong": "#2166ac", "Moderate": "#67a9cf", "Suggestive": "#d1e5f0"}
C_MIXED, C_NULL, C_NA, C_REFRAMES = "#9e6bb0", "#d9d9d9", "#fbfbfb", "#e6b800"

TOKEN_DISP = {
    "strong+": "Strong +", "mod+": "Mod +", "sugg+": "Sugg +",
    "strong-": "Strong −", "mod-": "Mod −", "sugg-": "Sugg −",
    "mixed": "mixed", "null": "null", "reframes": "reframes", "n.a.": "n.a.",
}


def _lum(hex_c):
    r, g, b = (int(hex_c[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.587 * b


def cell_style(cv, sign, grade):
    """(facecolor, textcolor, hatch, fontweight) for a cell."""
    if cv == "n.a.":
        return C_NA, "#9a9a9a", None, "normal"
    if cv == "null":
        return C_NULL, "#3a3a3a", None, "normal"
    if cv == "reframes":
        return C_REFRAMES, "#000000", "////", "bold"
    if cv == "mixed":
        return C_MIXED, "#ffffff", None, "normal"
    ramp = POS if sign == "+" else NEG
    fc = ramp.get(grade, ramp["Moderate"])
    tc = "#ffffff" if _lum(fc) < 0.5 else "#1a1a1a"
    return fc, tc, None, ("bold" if grade == "Strong" else "normal")


# ---- figure scaffold -------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "axes.linewidth": 0.8,
})
fig = plt.figure(figsize=(13.0, 8.6))
outer = fig.add_gridspec(1, 2, width_ratios=[1.5, 1.0], wspace=0.06,
                         left=0.205, right=0.985, top=0.79, bottom=0.07)
axA = fig.add_subplot(outer[0, 0])
right = outer[0, 1].subgridspec(2, 1, height_ratios=[1.05, 1.0], hspace=0.36)
axB = fig.add_subplot(right[0])
axL = fig.add_subplot(right[1])
axL.axis("off")

# ---- Panel A: convergence heatmap -----------------------------------------
GAP = 1.0                       # blank slot between mouse block and human band
human_slot = len(MOUSE_ROWS) + GAP
rows = [(name, arc, i) for i, (name, arc) in enumerate(MOUSE_ROWS)]
rows.append((HUMAN_ROW[0], HUMAN_ROW[1], human_slot))

# faint background band behind the human row to reinforce separation
axA.add_patch(Rectangle((-0.02, human_slot - 0.06), 3.04, 1.02,
                        facecolor="#f3f0f7", edgecolor="none", zorder=0))

for name, arc, slot in rows:
    for c, ax_key in enumerate(AXES):
        r = cell_of[(name, ax_key)]
        fc, tc, hatch, fw = cell_style(r.cell_value, r.sign, r.grade)
        axA.add_patch(Rectangle((c, slot), 1, 1, facecolor=fc, hatch=hatch,
                                edgecolor="white", linewidth=2.2, zorder=2))
        axA.text(c + 0.5, slot + 0.5, TOKEN_DISP.get(r.cell_value, r.cell_value),
                 ha="center", va="center", fontsize=8.3, color=tc,
                 fontweight=fw, zorder=3)

# divider + band caption between the mouse block and the human row
div_y = len(MOUSE_ROWS) + GAP / 2
axA.axhline(div_y, color="#888888", lw=1.0, ls=(0, (5, 3)), zorder=1)
axA.text(1.5, div_y - 0.16, "mouse: 11 orthogonal evidence layers + 4 base modalities",
         ha="center", va="bottom", fontsize=7.6, color="#666666", style="italic")
axA.text(1.5, div_y + 0.18,
         "human: separate band (directional translation; never counted in the contests)",
         ha="center", va="top", fontsize=7.6, color="#7a5aa0", style="italic")

axA.set_xlim(-0.03, 3.03)
axA.set_ylim(human_slot + 1.15, -0.12)         # inverted: row 0 at top
axA.set_xticks([0.5, 1.5, 2.5])
axA.set_xticklabels(AXIS_LABELS, fontsize=9.5, fontweight="bold")
axA.xaxis.set_ticks_position("top")
axA.tick_params(top=False, left=False)
ROW_DISPLAY = {"NF-kB attenuation": "NF-κB attenuation"}
yt = [slot + 0.5 for *_unused, slot in rows]
ylab = [f"{ROW_DISPLAY.get(name, name)}  ({arc})" for name, arc, _ in rows]
axA.set_yticks(yt)
axA.set_yticklabels(ylab, fontsize=8.7)
# italicise + recolour the human tick label
axA.get_yticklabels()[-1].set_style("italic")
axA.get_yticklabels()[-1].set_color("#7a5aa0")
for s in axA.spines.values():
    s.set_visible(False)

# ---- Panel B: contest margins ---------------------------------------------
cmap_axis = {
    "amyloid_activation": "Amyloid activation",
    "synaptic_suppression": "Synaptic suppression",
    "interaction_metabolic": "Tau×amyloid interaction",
}
con_sorted = con.sort_values("margin", ascending=True).reset_index(drop=True)
mmax = con_sorted["margin"].max()
ypos = list(range(len(con_sorted)))
axB.barh(ypos, con_sorted["margin"], height=0.5, color="#34648a",
         edgecolor="white", zorder=3)
axB.set_yticks([])
axB.set_ylim(-0.7, len(con_sorted) - 0.3)
axB.set_xlim(0, mmax * 1.16)
axB.set_xlabel("Contest margin  (winner − loser net-support)", fontsize=8.5)
for i, row in con_sorted.iterrows():
    win_id = row["winner"].split(" (")[0]
    win_name = row["winner"].split("(", 1)[1].rstrip(")")
    axB.text(0, i + 0.33, cmap_axis[row["contest"]], va="bottom", ha="left",
             fontsize=9.2, fontweight="bold", color="#222222")
    axB.text(row["margin"] + mmax * 0.012, i, f"margin {int(row['margin'])}",
             va="center", ha="left", fontsize=9.3, fontweight="bold",
             color="#222222")
    axB.text(0, i - 0.33,
             f"{win_id} wins · S/M/Su {int(row['n_strong'])}/{int(row['n_moderate'])}/{int(row['n_suggestive'])}"
             f" · {int(row['n_layers'])} layers - {win_name}",
             va="top", ha="left", fontsize=6.8, color="#555555")
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.spines["left"].set_visible(False)
axB.tick_params(left=False)

# ---- legend ----------------------------------------------------------------
legend_items = [
    Patch(facecolor=POS["Strong"], label="Strong  +"),
    Patch(facecolor=POS["Moderate"], label="Moderate  +"),
    Patch(facecolor=POS["Suggestive"], label="Suggestive  +"),
    Patch(facecolor=NEG["Strong"], label="Strong  −"),
    Patch(facecolor=NEG["Moderate"], label="Moderate  −"),
    Patch(facecolor=NEG["Suggestive"], label="Suggestive  −"),
    Patch(facecolor=C_MIXED, label="mixed  (sign-divergent)"),
    Patch(facecolor=C_NULL, label="null  (tested, no effect)"),
    Patch(facecolor=C_REFRAMES, hatch="////", label="reframes  (recasts the question)"),
    Patch(facecolor=C_NA, edgecolor="#cccccc", label="n.a.  (axis not probed)"),
]
leg = axL.legend(handles=legend_items, loc="upper left", ncol=2, frameon=False,
                 fontsize=8.2, handlelength=1.3, handleheight=1.2,
                 columnspacing=1.4, labelspacing=0.62,
                 title="Cell = sign × confidence grade   (token repeats grade for grayscale / colourblind reading)",
                 title_fontsize=8.4, bbox_to_anchor=(0.0, 1.02))
leg.get_title().set_fontweight("bold")
axL.text(0.0, -0.04,
         "Each cell cites a verdict-TSV value or ledger claim_id (capstone_convergence_matrix.tsv `source`).\n"
         "Amyloid activation is the strongly convergent axis (the lone − is tau's NF-κB attenuation, arc I).\n"
         "The interaction is mechanism-signed (TF Myc− / kinase Gsk3b+) and recovered in dynamics (M, O) &\n"
         "causal topology (J), yet null in composition (L), de-novo regulons (K) and human-under-FDR, and is\n"
         "NOT microglia-specific (N): convergence with honest, first-class nulls.",
         ha="left", va="top", fontsize=7.5, color="#333333", transform=axL.transAxes)

fig.suptitle("Integrated convergence of the tau×amyloid microglial model across 12 evidence layers + human translation",
             x=0.205, ha="left", fontsize=12.4, fontweight="bold", y=0.965)
fig.text(0.205, 0.928,
         "Dual view: (A) each arc's own signed verdict per biological axis (the landscape: orthogonal corroboration the contest arithmetic discounts);  "
         "(B) the ledger-adjudicated margins that moved the verdicts.",
         ha="left", fontsize=8.6, color="#444444")
# panel labels on one baseline, above the column-header zone
fig.text(0.205, 0.852, "A   Landscape view: each arc's signed verdict per axis",
         ha="left", fontsize=10.8, fontweight="bold")
fig.text(axB.get_position().x0, 0.852,
         "B   Contest view: adjudicated margins (sealed §17)",
         ha="left", fontsize=10.8, fontweight="bold")

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=300, facecolor="white", bbox_inches="tight")
print(f"[capstone-fig] wrote {OUT}  ({fig.get_size_inches()[0]:.1f}x{fig.get_size_inches()[1]:.1f} in @300dpi)")
