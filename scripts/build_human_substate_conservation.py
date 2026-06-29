#!/usr/bin/env python3
"""scripts/build_human_substate_conservation.py

Session H3 of the human cross-species validation plan
(storage/notes/human_validation_plan.md).

Question (H3): do the four MOUSE microglial substates -- homeostatic / DAM /
IFN / proliferative -- resolve as distinct populations in human SEA-AD
microglia, and how do they map onto the dataset-native SEA-AD `Supertype`
labels? This gates whether H4 may make per-state interaction claims
(anti-anchoring #3: conservation precedes per-state claims).

Method (faithful reproduction of the mouse labelling):
  The mouse pipeline labelled every nucleus by scoring the four canonical
  marker lists (R/constants.R::canonical_microglia_markers) with
  Seurat::AddModuleScore(nbin = 12, ctrl = 50) on the SCT assay and taking
  the per-nucleus argmax (R/microglia.R::label_microglia_states). There is
  NO de-novo per-state marker table; these curated lists ARE the substate
  definitions. We therefore score the SAME four lists (mapped to human
  orthologs in H2) per human cell with scanpy's score_genes -- the direct
  AddModuleScore analogue -- matching n_bins = 12, ctrl_size = 50, and take
  the argmax. score_genes is run on X (CELLxGENE log1p-normalised) with
  use_raw=False; .raw holds raw COUNTS and would be wrong to score on.
  Residual nuance: the human X is library-size log1p, not SCT (we do not
  re-SCT the human data); this is the standard cross-dataset scoring setup.

  All 26 H2 signatures are scored in the same pass (the h5ad load is the
  expensive step), so H4/H5 inherit the per-cell + per-donon-x-state score
  panel without re-loading the data. This extends H3's stated 4-signature
  scope for efficiency; the conservation analysis itself uses only the four
  substate signatures.

Conservation readout: argmax predicted_substate x SEA-AD Supertype cross-
tab, with adjusted Rand index, (adjusted) normalised mutual information,
Cramer's V, and a per-state resolution verdict (cell fraction, mean self-
score, dominant Supertype enrichment). The verdict for DAM and IFN is the
H4 gate.

Inputs:
  storage/data/seaad/microglia_{mtg,dlpfc}.h5ad        (H1)
  storage/cache/human_validation_signatures.rds        (H2; via JSON bridge)

Outputs (storage/cache/ = heavier, H4/H5-consumed):
  human_validation_signatures_human.json   RDS->JSON bridge (human lists + meta)
  human_substate_percell.csv.gz            per cell: 26 scores, labels, UMAP
  human_substate_score_means.csv           per donor x region x state: mean scores (primary H4 input)
  human_substate_pseudobulk_counts.csv.gz  genes x sample: raw-count sums (signature-union genes)
  human_substate_pseudobulk_samples.csv    sample metadata for the pseudobulk matrix
  human_substate_conservation_provenance.txt
Outputs (storage/results/ = human-readable summaries):
  human_substate_composition.tsv           per donor x region: predicted + supertype fractions
  human_substate_crosstab.tsv              predicted_substate x Supertype counts (overall + per region)
  human_substate_conservation_metrics.tsv  ARI / AMI / NMI / Cramer's V + per-state verdict
  human_substate_signature_coverage.tsv    data-universe intersection (n signature genes present in SEA-AD var)

Run: ./.venv/bin/python scripts/build_human_substate_conservation.py [--smoke N]
"""
import os
import sys
import json
import subprocess
import textwrap
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             adjusted_mutual_info_score)

sc.settings.verbosity = 0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEAAD = os.path.join(ROOT, "storage", "data", "seaad")
CACHE = os.path.join(ROOT, "storage", "cache")
RESULTS = os.path.join(ROOT, "storage", "results")

SOURCES = {
    "MTG":   os.path.join(SEAAD, "microglia_mtg.h5ad"),
    "DLPFC": os.path.join(SEAAD, "microglia_dlpfc.h5ad"),
}

# substate -> its H2 signature name; order defines the argmax label vector.
SUBSTATES = ["homeostatic", "DAM", "IFN", "proliferative"]
SUBSTATE_SIG = {s: f"substate_{s}" for s in SUBSTATES}

# match Seurat::AddModuleScore(nbin = 12, ctrl = 50) from the mouse pipeline.
N_BINS, CTRL, RANDOM_STATE = 12, 50, 0

SIG_RDS = os.path.join(CACHE, "human_validation_signatures.rds")
SIG_JSON = os.path.join(CACHE, "human_validation_signatures_human.json")


def load_signatures():
    """Bridge the H2 RDS into Python via a jsonlite dump (pyreadr absent)."""
    r_code = (
        f'suppressWarnings(suppressMessages(library(jsonlite)));'
        f's <- readRDS("{SIG_RDS}");'
        f'jsonlite::write_json(list(human = s$human, meta = s$meta),'
        f' "{SIG_JSON}", auto_unbox = FALSE, dataframe = "rows", pretty = TRUE)'
    )
    subprocess.run(["Rscript", "-e", r_code], check=True,
                   cwd=ROOT, stdout=subprocess.DEVNULL)
    with open(SIG_JSON) as fh:
        obj = json.load(fh)
    human = {k: [str(g) for g in (v if isinstance(v, list) else [v])]
             for k, v in obj["human"].items()}
    meta = obj["meta"]
    return human, meta


def score_region(region, path, human, smoke=0):
    """Load one region, score all signatures, argmax substate; return
    (per-cell DataFrame, raw-count submatrix DataFrame over union genes,
     coverage rows)."""
    a = sc.read_h5ad(path)
    if smoke:
        a = a[:smoke].copy()
    # symbol space: var_names are Ensembl, symbols in feature_name.
    a.var["ensembl"] = a.var_names.astype(str)
    a.var_names = a.var["feature_name"].astype(str)
    a.var_names_make_unique()
    var_set = set(a.var_names)

    sig_names = list(human.keys())
    coverage = []
    for s in sig_names:
        present = [g for g in human[s] if g in var_set]
        coverage.append(dict(signature=s, region=region,
                             n_human=len(human[s]), n_present=len(present),
                             pct_present=round(100 * len(present) /
                                               max(len(human[s]), 1), 1)))
        col = f"score_{s}"
        if not present:
            a.obs[col] = np.nan
            continue
        sc.tl.score_genes(a, present, score_name=col, ctrl_size=CTRL,
                          n_bins=N_BINS, random_state=RANDOM_STATE,
                          use_raw=False)

    # argmax over the four substate scores -> predicted_substate
    sub_cols = [f"score_{SUBSTATE_SIG[s]}" for s in SUBSTATES]
    ms = a.obs[sub_cols].to_numpy(dtype=float)
    pred = np.array(SUBSTATES)[np.argmax(ms, axis=1)]

    umap = a.obsm["X_umap"]
    df = pd.DataFrame({
        "cell_id": a.obs_names.astype(str),
        "donor_id": a.obs["donor_id"].astype(str).to_numpy(),
        "region": region,
        "is_reference": a.obs["Thal phase"].astype(str).eq("Reference").to_numpy(),
        "Supertype": a.obs["Supertype"].astype(str).to_numpy(),
        "Subclass": a.obs["Subclass"].astype(str).to_numpy(),
        "predicted_substate": pred,
        "UMAP1": umap[:, 0], "UMAP2": umap[:, 1],
    })
    for s in sig_names:
        df[f"score_{s}"] = a.obs[f"score_{s}"].to_numpy()

    # ---- raw-count pseudobulk submatrix over signature-union genes --------
    union_genes = sorted(set(g for v in human.values() for g in v))
    present_union = [g for g in union_genes if g in var_set]
    raw = a.raw.to_adata()
    # raw and main share gene order (same 35,483-gene space); re-label raw to
    # the made-unique symbols so a symbol subset is consistent with scoring.
    assert raw.n_vars == a.n_vars, "raw/main gene-count mismatch"
    assert list(raw.var_names) == list(a.var["ensembl"]), \
        "raw/main gene-order mismatch (cannot positionally re-label)"
    raw.var_names = a.var_names
    sub = raw[:, present_union]
    M = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    expr = pd.DataFrame(M, columns=present_union)
    expr.insert(0, "predicted_substate", pred)
    expr.insert(0, "donor_id", df["donor_id"].to_numpy())
    return df, expr, coverage


def pseudobulk_counts(expr_all, gene_cols):
    """Sum raw counts per (donor, region, state) and per (donor, region, all);
    return wide genes x sample matrix + sample metadata."""
    mat_cols, sample_rows = {}, []

    def add(group_df, donor, region, state):
        key = f"{donor}|{region}|{state}"
        mat_cols[key] = group_df[gene_cols].sum(axis=0).to_numpy()
        sample_rows.append(dict(sample_key=key, donor_id=donor, region=region,
                                state=state, n_cells=len(group_df)))

    for (donor, region), g in expr_all.groupby(["donor_id", "region"],
                                               observed=True, sort=True):
        add(g, donor, region, "all")
        for state, gs in g.groupby("predicted_substate", observed=True, sort=True):
            if len(gs):
                add(gs, donor, region, state)

    mat = pd.DataFrame(mat_cols, index=gene_cols)
    mat.index.name = "gene"
    samples = pd.DataFrame(sample_rows)
    return mat, samples


def cramers_v(ct):
    """Cramer's V from a contingency-count DataFrame."""
    obs = ct.to_numpy(dtype=float)
    n = obs.sum()
    if n == 0:
        return np.nan
    row, col = obs.sum(1, keepdims=True), obs.sum(0, keepdims=True)
    exp = row @ col / n
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.nansum(np.where(exp > 0, (obs - exp) ** 2 / exp, 0.0))
    k = min(obs.shape) - 1
    return float(np.sqrt(chi2 / (n * k))) if k > 0 else np.nan


def main():
    smoke = 0
    if "--smoke" in sys.argv:
        i = sys.argv.index("--smoke")
        smoke = int(sys.argv[i + 1]) if i + 1 < len(sys.argv) else 2000

    human, meta = load_signatures()
    sig_names = list(human.keys())
    score_cols = [f"score_{s}" for s in sig_names]
    union_genes = sorted(set(g for v in human.values() for g in v))
    print(f"[H3] {len(sig_names)} signatures | "
          f"{len(union_genes)} union genes | smoke={smoke or 'OFF'}")

    percell, exprs, coverage = [], [], []
    for region, path in SOURCES.items():
        if not os.path.exists(path):
            sys.exit(f"MISSING input h5ad: {path}")
        df, expr, cov = score_region(region, path, human, smoke=smoke)
        percell.append(df)
        exprs.append(expr)
        coverage += cov
        print(f"  [{region}] {len(df)} cells scored | "
              f"predicted: {dict(df['predicted_substate'].value_counts())}")

    pc = pd.concat(percell, ignore_index=True)
    expr_all = pd.concat(exprs, ignore_index=True)
    expr_all["region"] = pc["region"].to_numpy()
    gene_cols = [g for g in union_genes if g in expr_all.columns]

    cov = pd.DataFrame(coverage)

    # ---- conservation: cross-tab + metrics --------------------------------
    def metrics_block(sub_pc, label):
        ari = adjusted_rand_score(sub_pc["Supertype"], sub_pc["predicted_substate"])
        ami = adjusted_mutual_info_score(sub_pc["Supertype"], sub_pc["predicted_substate"])
        nmi = normalized_mutual_info_score(sub_pc["Supertype"], sub_pc["predicted_substate"])
        ct = pd.crosstab(sub_pc["predicted_substate"], sub_pc["Supertype"])
        return dict(scope=label, n_cells=len(sub_pc), ARI=round(ari, 4),
                    AMI=round(ami, 4), NMI=round(nmi, 4),
                    cramers_v=round(cramers_v(ct), 4))

    metric_rows = [metrics_block(pc, "all_regions")]
    for region in SOURCES:
        metric_rows.append(metrics_block(pc[pc.region == region], region))

    # cross-tab (long) overall + per region
    ct_rows = []
    for scope, sub in [("all_regions", pc)] + [(r, pc[pc.region == r]) for r in SOURCES]:
        ct = pd.crosstab(sub["predicted_substate"], sub["Supertype"])
        for ps in ct.index:
            for st in ct.columns:
                ct_rows.append(dict(scope=scope, predicted_substate=ps,
                                    Supertype=st, n=int(ct.loc[ps, st])))
    crosstab = pd.DataFrame(ct_rows)

    # ---- per-state resolution verdict (the H4 gate) -----------------------
    # SEA-AD splits microglia into a dominant homeostatic-leaning supertype
    # (Micro-PVM_2, ~75-85% of cells) plus rarer disease-emergent "-SEAAD"
    # supertypes. "Most common supertype" is therefore Micro-PVM_2 for every
    # predicted state and uninformative; the conservation signal is whether a
    # state ENRICHES in the disease "-SEAAD" supertypes (fold over global) and
    # what fraction of its cells are in any "-SEAAD" supertype.
    global_super = pc["Supertype"].value_counts(normalize=True)
    is_seaad = pc["Supertype"].str.endswith("-SEAAD")
    global_seaad = float(is_seaad.mean())
    n_total = len(pc)
    MIN_CT = 20  # min cells in a state x supertype cell to trust a fold
    verdict_rows = []
    for s in SUBSTATES:
        cells = pc[pc["predicted_substate"] == s]
        frac = len(cells) / n_total
        mean_self = float(cells[f"score_{SUBSTATE_SIG[s]}"].mean()) if len(cells) else np.nan
        frac_seaad = float(cells["Supertype"].str.endswith("-SEAAD").mean()) if len(cells) else np.nan
        # most-enriched supertype: max fold over its global share, guarded by count
        enr_super, enr_fold = "", np.nan
        if len(cells):
            shares = cells["Supertype"].value_counts(normalize=True)
            counts = cells["Supertype"].value_counts()
            folds = {st: shares[st] / global_super.get(st, np.nan)
                     for st in shares.index if counts[st] >= MIN_CT}
            if folds:
                enr_super = max(folds, key=folds.get)
                enr_fold = float(folds[enr_super])
        resolves = bool(frac >= 0.005 and mean_self > 0)
        verdict_rows.append(dict(
            predicted_substate=s, n_cells=len(cells), frac=round(frac, 4),
            mean_self_score=round(mean_self, 4),
            frac_SEAAD_supertype=round(frac_seaad, 4),
            SEAAD_fold_vs_global=round(frac_seaad / global_seaad, 3)
            if global_seaad > 0 else np.nan,
            most_enriched_supertype=enr_super,
            enrichment_fold=round(enr_fold, 3) if enr_fold == enr_fold else np.nan,
            resolves=resolves))
    verdict = pd.DataFrame(verdict_rows)

    # ---- per-donor x region x state score means (primary H4 input) --------
    grp = pc.groupby(["donor_id", "region", "predicted_substate"], observed=True)
    means_state = grp[score_cols].mean()
    means_state["n_cells"] = grp.size()
    means_state = means_state.reset_index().rename(
        columns={"predicted_substate": "state"})
    grp_all = pc.groupby(["donor_id", "region"], observed=True)
    means_all = grp_all[score_cols].mean()
    means_all["n_cells"] = grp_all.size()
    means_all = means_all.reset_index()
    means_all["state"] = "all"
    score_means = pd.concat([means_all, means_state], ignore_index=True)
    score_means = score_means[["donor_id", "region", "state", "n_cells"] + score_cols]

    # ---- per-donor x region composition (predicted + supertype) -----------
    comp_rows = []
    for (donor, region), g in pc.groupby(["donor_id", "region"], observed=True):
        n = len(g)
        for lab, vc in [("predicted", g["predicted_substate"].value_counts()),
                        ("supertype", g["Supertype"].value_counts())]:
            for label, cnt in vc.items():
                comp_rows.append(dict(donor_id=donor, region=region,
                                      label_type=lab, label=label,
                                      n_cells=int(cnt), frac=round(cnt / n, 4)))
    composition = pd.DataFrame(comp_rows)

    # ---- pseudobulk count matrix ------------------------------------------
    pb_mat, pb_samples = pseudobulk_counts(expr_all, gene_cols)

    # ====================== write outputs ==================================
    out = {}
    if smoke:
        print("\n[SMOKE] skipping writes; shape summary only:")
        print("  per-cell:", pc.shape, "| score_means:", score_means.shape,
              "| pb_mat:", pb_mat.shape, "| pb_samples:", pb_samples.shape)
        print("\nconservation metrics:")
        print(pd.DataFrame(metric_rows).to_string(index=False))
        print("\nper-state verdict:")
        print(verdict.to_string(index=False))
        return

    pc.to_csv(os.path.join(CACHE, "human_substate_percell.csv.gz"),
              index=False, compression="gzip")
    score_means.to_csv(os.path.join(CACHE, "human_substate_score_means.csv"),
                       index=False)
    pb_mat.to_csv(os.path.join(CACHE, "human_substate_pseudobulk_counts.csv.gz"),
                  compression="gzip")
    pb_samples.to_csv(os.path.join(CACHE, "human_substate_pseudobulk_samples.csv"),
                      index=False)
    composition.to_csv(os.path.join(RESULTS, "human_substate_composition.tsv"),
                       sep="\t", index=False)
    crosstab.to_csv(os.path.join(RESULTS, "human_substate_crosstab.tsv"),
                    sep="\t", index=False)
    cov.to_csv(os.path.join(RESULTS, "human_substate_signature_coverage.tsv"),
               sep="\t", index=False)
    metrics = pd.DataFrame(metric_rows)
    with open(os.path.join(RESULTS, "human_substate_conservation_metrics.tsv"), "w") as fh:
        fh.write("# clustering-agreement metrics (predicted_substate vs SEA-AD Supertype)\n")
        metrics.to_csv(fh, sep="\t", index=False)
        fh.write("\n# per-state resolution verdict (DAM & IFN rows are the H4 gate)\n")
        verdict.to_csv(fh, sep="\t", index=False)

    # ---- provenance + gate narrative --------------------------------------
    low_cov = cov[cov.pct_present < 60].sort_values("pct_present")
    dam_ok = bool(verdict.set_index("predicted_substate").loc["DAM", "resolves"])
    ifn_ok = bool(verdict.set_index("predicted_substate").loc["IFN", "resolves"])
    lines = []
    lines.append("SEA-AD human microglial substate conservation -- H3")
    lines.append("=" * 72)
    lines.append(f"Cells scored: {len(pc)} ({', '.join(f'{r}={int((pc.region==r).sum())}' for r in SOURCES)})")
    lines.append(f"Donors: {pc['donor_id'].nunique()} | reference (neurotypical) cells: {int(pc.is_reference.sum())}")
    lines.append("")
    lines.append("Scoring: scanpy score_genes on X (log1p), use_raw=False, "
                 f"n_bins={N_BINS}, ctrl_size={CTRL}, random_state={RANDOM_STATE} "
                 "-- matches mouse Seurat::AddModuleScore(nbin=12, ctrl=50); "
                 "predicted_substate = argmax of the four substate scores "
                 "(reproduces R/microglia.R::label_microglia_states).")
    lines.append("")
    lines.append("Clustering agreement (predicted_substate vs SEA-AD Supertype):")
    lines.append(textwrap.indent(metrics.to_string(index=False), "  "))
    lines.append("")
    lines.append(f"Global '-SEAAD' (disease-emergent) supertype fraction: {global_seaad:.4f}")
    lines.append("Per-state resolution verdict (cell fraction; mean self-score; "
                 "fraction in any '-SEAAD' supertype + fold over global; "
                 "most-enriched supertype + fold, min 20 cells):")
    lines.append(textwrap.indent(verdict.to_string(index=False), "  "))
    lines.append("")
    lines.append("Cross-tab (all regions): predicted_substate x Supertype")
    lines.append(textwrap.indent(
        pd.crosstab(pc["predicted_substate"], pc["Supertype"]).to_string(), "  "))
    lines.append("")
    lines.append("GATE (H4 per-state claims, anti-anchoring #3):")
    lines.append(f"  DAM resolves as a distinct human population: {dam_ok}")
    lines.append(f"  IFN resolves as a distinct human population: {ifn_ok}")
    lines.append("  Rule: resolves := (cell fraction >= 0.5% of microglia) AND "
                 "(mean self-score > 0). Supertype enrichment is reported as "
                 "corroborating biology, not part of the pass/fail rule. The "
                 "verdict is advisory for H4: a state that fails is restricted "
                 "to whole-microglia modelling; a null per-state interaction on "
                 "a resolving state remains interpretable under collinearity.")
    lines.append("")
    lines.append("Signature data-universe coverage (genes present in SEA-AD var; "
                 "the H2-deferred intersection). Full table: "
                 "storage/results/human_substate_signature_coverage.tsv.")
    if len(low_cov):
        lines.append(f"  {len(low_cov)} signature(s) < 60% present in SEA-AD var:")
        for _, r in low_cov.iterrows():
            lines.append(f"    {r.signature} [{r.region}]: "
                         f"{int(r.n_present)}/{int(r.n_human)} ({r.pct_present}%)")
    else:
        lines.append("  all signatures >= 60% present in SEA-AD var.")
    lines.append("")
    lines.append("Outputs:")
    lines.append("  cache/human_substate_percell.csv.gz           per-cell 26 scores + labels + UMAP")
    lines.append("  cache/human_substate_score_means.csv          per donor x region x state mean scores (H4 primary)")
    lines.append("  cache/human_substate_pseudobulk_counts.csv.gz genes x sample raw-count sums (union genes)")
    lines.append("  cache/human_substate_pseudobulk_samples.csv   pseudobulk sample metadata")
    lines.append("  results/human_substate_composition.tsv        per donor x region predicted+supertype fractions")
    lines.append("  results/human_substate_crosstab.tsv           predicted x Supertype counts (overall+region)")
    lines.append("  results/human_substate_conservation_metrics.tsv  ARI/AMI/NMI/Cramer's V + per-state verdict")
    lines.append("  results/human_substate_signature_coverage.tsv  signature gene presence in SEA-AD var")
    with open(os.path.join(CACHE, "human_substate_conservation_provenance.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- smoke summary to stdout ------------------------------------------
    print("\n=== conservation metrics ===")
    print(metrics.to_string(index=False))
    print("\n=== per-state verdict (DAM & IFN = H4 gate) ===")
    print(verdict.to_string(index=False))
    print(f"\nGATE: DAM resolves={dam_ok} | IFN resolves={ifn_ok}")
    print(f"\nWrote per-cell ({pc.shape}), score_means ({score_means.shape}), "
          f"pseudobulk ({pb_mat.shape} genes x samples), + 5 result TSVs.")


if __name__ == "__main__":
    main()
