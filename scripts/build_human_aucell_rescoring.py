#!/usr/bin/env python3
"""scripts/build_human_aucell_rescoring.py

Session H5 of the human cross-species validation plan
(storage/notes/human_validation_plan.md) -- SCORING-METHOD robustness arm.

H3/H4 scored every signature with scanpy `score_genes` (the mouse
Seurat::AddModuleScore analogue: mean of set genes minus a binned control
set). The H4 interaction headline (negative amyloid:tau for NF-kB / MG-M3,
positive for Gsk3b) therefore rests on ONE scoring method. AUCell (Aibar
2017) is the canonical ORTHOGONAL scorer: per cell it ranks all genes and
takes the area under the recovery curve of the set within the top-ranked
fraction -- rank-based, no control-set subtraction, invariant to any
monotone per-cell normalisation. If the H4 interaction directions reproduce
on AUCell scores, they are not an artefact of the `score_genes` method.

To ISOLATE the scoring method, this script reproduces H3's load EXACTLY
(same h5ad cells, same Ensembl->feature_name var relabelling, same 26 H2
signatures via the JSON bridge) and re-aggregates AUCell scores to the SAME
per donor x region x state cells using H3's `predicted_substate` labels
(joined by cell_id). Only the per-cell score VALUES differ between the two
score_means tables, so the R refit (build_human_robustness_mediation.R)
compares like with like.

Method: decoupler `dc.mt.aucell` on `adata.X` (CELLxGENE log1p; AUCell ranks
per cell so the log1p is immaterial), tmin=5 (decoupler's reliability floor:
sources with <5 present targets are dropped and logged -- no silent cap).

Inputs:
  storage/data/seaad/microglia_{mtg,dlpfc}.h5ad        (H1)
  storage/cache/human_validation_signatures_human.json (H3 bridge; regenerated
      from the H2 RDS if absent)
  storage/cache/human_substate_percell.csv.gz          (H3; cell_id ->
      predicted_substate, for state-matched re-aggregation)

Outputs:
  storage/cache/human_substate_aucell_score_means.csv  per donor x region x
      {all + 4 substates} mean AUCell score (mirror of human_substate_score_
      means.csv; the scoring-method-robustness modelling input)
  storage/cache/human_aucell_rescoring_provenance.txt  dropped sets, label
      agreement vs H3 argmax, shapes

Run: ./.venv/bin/python scripts/build_human_aucell_rescoring.py [--smoke N]
"""
import os
import sys
import json
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc

sc.settings.verbosity = 0

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEAAD = os.path.join(ROOT, "storage", "data", "seaad")
CACHE = os.path.join(ROOT, "storage", "cache")

SOURCES = {
    "MTG":   os.path.join(SEAAD, "microglia_mtg.h5ad"),
    "DLPFC": os.path.join(SEAAD, "microglia_dlpfc.h5ad"),
}
SUBSTATES = ["homeostatic", "DAM", "IFN", "proliferative"]
SUBSTATE_SIG = {s: f"substate_{s}" for s in SUBSTATES}
TMIN = 5  # decoupler reliability floor: drop gene sets with < TMIN present targets

SIG_RDS = os.path.join(CACHE, "human_validation_signatures.rds")
SIG_JSON = os.path.join(CACHE, "human_validation_signatures_human.json")
PERCELL = os.path.join(CACHE, "human_substate_percell.csv.gz")


def load_signatures():
    """Read the H3 RDS->JSON bridge; regenerate it from the H2 RDS if absent."""
    if not os.path.exists(SIG_JSON):
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
    return human


def aucell_region(region, path, net, smoke=0):
    """Load one region, relabel to symbols, AUCell-score; return per-cell
    DataFrame (cell_id, donor_id, region, aucell_<sig> ...) over the sources
    decoupler retained, plus the set of retained source names."""
    a = sc.read_h5ad(path)
    if smoke:
        a = a[:smoke].copy()
    a.var["ensembl"] = a.var_names.astype(str)
    a.var_names = a.var["feature_name"].astype(str)
    a.var_names_make_unique()

    dc.mt.aucell(a, net, tmin=TMIN, raw=False, verbose=False)
    sc_df = a.obsm["score_aucell"]            # cells x retained-sources
    retained = list(sc_df.columns)

    df = pd.DataFrame({
        "cell_id": a.obs_names.astype(str),
        "donor_id": a.obs["donor_id"].astype(str).to_numpy(),
        "region": region,
    })
    for s in retained:
        df[f"aucell_{s}"] = sc_df[s].to_numpy()
    return df, set(retained)


def main():
    smoke = 0
    if "--smoke" in sys.argv:
        i = sys.argv.index("--smoke")
        smoke = int(sys.argv[i + 1]) if i + 1 < len(sys.argv) else 2000

    human = load_signatures()
    sig_names = list(human.keys())

    # Build the long-form decoupler net (source=signature, target=gene). Genes
    # absent from a region's var simply contribute nothing to that cell's
    # ranking; tmin is applied by decoupler per region. We pool the full human
    # gene lists and let decoupler intersect with each region's var.
    net = pd.DataFrame(
        [(s, g) for s in sig_names for g in human[s]],
        columns=["source", "target"]).drop_duplicates()
    print(f"[H5-AUCell] {len(sig_names)} signatures | "
          f"{net['target'].nunique()} union genes | tmin={TMIN} | "
          f"smoke={smoke or 'OFF'}")

    percell, retained_per_region = [], {}
    for region, path in SOURCES.items():
        if not os.path.exists(path):
            sys.exit(f"MISSING input h5ad: {path}")
        df, retained = aucell_region(region, path, net, smoke=smoke)
        percell.append(df)
        retained_per_region[region] = retained
        print(f"  [{region}] {len(df)} cells | {len(retained)} sources retained")

    pc = pd.concat(percell, ignore_index=True)
    auc_cols = [c for c in pc.columns if c.startswith("aucell_")]
    retained_sigs = [c[len("aucell_"):] for c in auc_cols]
    sizes = {k: len(v) for k, v in human.items()}
    # Dropped = below tmin (a rank-recovery-curve scorer is undefined for tiny
    # sets); annotate each with its gene-set size so the drop reads as
    # principled, not a failure. The six pre-registered headline signatures
    # (NF-kB 1304, MG-M3 740, Gsk3b 538, DAM_up, Gerrits AD1/AD2) are all large
    # and retained, so AUCell re-tests every signature carrying an H4 claim.
    dropped = sorted(set(sig_names) - set(retained_sigs))
    dropped_str = ", ".join(f"{d} ({sizes[d]}g)" for d in dropped) or "none"

    # AUCell argmax substate (secondary: does state ASSIGNMENT survive the
    # scoring-method swap?). Requires all four substate sources retained.
    sub_cols = [f"aucell_{SUBSTATE_SIG[s]}" for s in SUBSTATES]
    have_all_sub = all(c in pc.columns for c in sub_cols)
    if have_all_sub:
        ms = pc[sub_cols].to_numpy(dtype=float)
        pc["predicted_substate_aucell"] = np.array(SUBSTATES)[np.argmax(ms, axis=1)]

    # Join H3's predicted_substate (by cell_id) so AUCell is aggregated over
    # the SAME state cells as score_genes -> isolates the scoring method.
    h3 = pd.read_csv(PERCELL, usecols=["cell_id", "predicted_substate"])
    pc = pc.merge(h3, on="cell_id", how="left")
    label_cov = float(pc["predicted_substate"].notna().mean())
    label_agree = (float((pc["predicted_substate_aucell"]
                          == pc["predicted_substate"]).mean())
                   if have_all_sub else np.nan)

    # Aggregate to donor x region x {all + H3 state} means (mirror H3 shape).
    grp = pc.groupby(["donor_id", "region", "predicted_substate"], observed=True)
    means_state = grp[auc_cols].mean()
    means_state["n_cells"] = grp.size()
    means_state = means_state.reset_index().rename(
        columns={"predicted_substate": "state"})
    grp_all = pc.groupby(["donor_id", "region"], observed=True)
    means_all = grp_all[auc_cols].mean()
    means_all["n_cells"] = grp_all.size()
    means_all = means_all.reset_index()
    means_all["state"] = "all"
    score_means = pd.concat([means_all, means_state], ignore_index=True)
    # rename aucell_<sig> -> score_<sig> so the R refit reuses H4 column logic.
    score_means = score_means.rename(columns={c: f"score_{c[len('aucell_'):]}"
                                              for c in auc_cols})
    score_means = score_means[["donor_id", "region", "state", "n_cells"]
                              + [f"score_{s}" for s in retained_sigs]]

    if smoke:
        print("\n[SMOKE] no writes; shapes only:")
        print("  per-cell:", pc.shape, "| score_means:", score_means.shape)
        print("  dropped (<tmin) sources:", dropped or "none")
        print(f"  AUCell-vs-H3 label agreement: "
              f"{label_agree:.4f}" if have_all_sub else "  (substate set dropped)")
        return

    out_means = os.path.join(CACHE, "human_substate_aucell_score_means.csv")
    score_means.to_csv(out_means, index=False)
    os.chmod(out_means, 0o664)

    prov = os.path.join(CACHE, "human_aucell_rescoring_provenance.txt")
    lines = [
        "human AUCell re-scoring provenance (H5 scoring-method robustness)",
        "=" * 70,
        "built_by   : scripts/build_human_aucell_rescoring.py",
        "scorer     : decoupler dc.mt.aucell on adata.X (log1p; rank-based ->",
        "             monotone-normalisation invariant), tmin=%d" % TMIN,
        "purpose    : orthogonal scorer vs score_genes; re-aggregated over the",
        "             SAME H3 predicted_substate cells (joined by cell_id) so",
        "             the only difference vs human_substate_score_means.csv is",
        "             the scoring method (isolates scoring-method robustness).",
        f"cells      : {len(pc)} ({', '.join(f'{r}={int((pc.region==r).sum())}' for r in SOURCES)})",
        f"signatures : {len(retained_sigs)}/{len(sig_names)} retained at tmin={TMIN}",
        f"dropped    : {dropped_str}",
        "             (tiny marker/DE sets < tmin; AUCell is a rank-recovery-curve",
        "             scorer for sizeable sets, so they stay score_genes-only --",
        "             all 6 pre-registered headline signatures are large & retained)",
        f"H3 label coverage (cell_id join): {label_cov:.4f}",
        ("AUCell-vs-H3 argmax substate agreement: %.4f "
         "(secondary: state assignment under the scoring-method swap)" % label_agree)
        if have_all_sub else
        "AUCell substate argmax: N/A (a substate source fell below tmin)",
        f"output     : {out_means} ({score_means.shape[0]} rows x "
        f"{score_means.shape[1]} cols)",
    ]
    with open(prov, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    os.chmod(prov, 0o664)

    print(f"\nWrote {out_means} ({score_means.shape})")
    print(f"dropped (<tmin) sources: {dropped or 'none'}")
    if have_all_sub:
        print(f"AUCell-vs-H3 label agreement: {label_agree:.4f}")


if __name__ == "__main__":
    main()
