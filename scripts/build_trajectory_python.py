#!/usr/bin/env python
# --------------------------------------------------------------------
# Trajectory arc (plan M3) -- the PYTHON half of the maximal-triangulation
# pseudotime build. Slingshot (R, primary) is computed by build_trajectory.R;
# this script adds the two INDEPENDENT cross-method partners the concordance
# check (guardrail #6) needs, both VELOCITY-FREE (no spliced layers exist):
#   * scanpy PAGA + DPT  -> a diffusion pseudotime that does NOT use Slingshot's
#     principal curves, so Slingshot-vs-DPT Spearman is a genuine method check.
#   * CellRank2 (cellrank 2.x) PseudotimeKernel -> GPCCA -> DAM-fate absorption
#     probability: re-expresses the readout as a fate probability (the M1 gate's
#     option-2 flavour) rather than a 1-D ordering. Best-effort / non-blocking.
#
# It is a thin, language-neutral bridge mirroring the SCENIC R->Python export
# pattern (build_scenic_grn.py): R owns the Seurat .rds + harmony embedding +
# the runtime-derived `state` labels and writes a compact CSV bundle; this owns
# scanpy + cellrank and writes a per-cell CSV back. NO expression matrix is
# needed -- PAGA/DPT/CellRank all run on the harmony embedding + the 4 states,
# the SAME low-dim space Slingshot uses, so the three methods differ only in
# algorithm, not in input (a clean concordance design).
#
# Inputs  (export_dir, written by build_trajectory.R):
#   embedding.csv   index = cell_id, columns = harmony_1..harmony_<n_dims>
#   obs.csv         cell_id, state, genotype, batch, genotype_batch
# Outputs (export_dir):
#   python_pseudotime.csv   cell_id, dpt_pseudotime, cellrank_dam_fate
#   python_provenance.json  root cell, PAGA connectivities, CellRank status,
#                           versions (human + LLM readable)
#
# Run (project .venv): .venv/bin/python scripts/build_trajectory_python.py \
#          --export-dir storage/cache/trajectory
# --------------------------------------------------------------------
import argparse
import json
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.simplefilter("ignore")


def log(msg):
    print(f"[build_trajectory_python] {msg}", flush=True)


def pick_root_cell(adata, root_state, end_state, dc_key="X_diffmap"):
    """Deterministic root = the `root_state` cell furthest (in diffusion space)
    from the `end_state` centroid -- the manifold tip of the root cluster, a
    principled DPT root (vs an arbitrary first cell). Falls back to the first
    root_state cell if the end_state is absent."""
    states = adata.obs["state"].to_numpy()
    dc = adata.obsm[dc_key]
    root_idx = np.flatnonzero(states == root_state)
    if root_idx.size == 0:
        raise ValueError(f"no '{root_state}' cells to root on")
    end_idx = np.flatnonzero(states == end_state)
    if end_idx.size == 0:
        return int(root_idx[0])
    end_centroid = dc[end_idx].mean(axis=0)
    d = np.linalg.norm(dc[root_idx] - end_centroid, axis=1)
    return int(root_idx[np.argmax(d)])


def run_cellrank(adata, dam_state, side_state, seed, auto=False):
    """Velocity-free CellRank2 DAM-fate probability via the CHEAP, PRINCIPLED
    route. The biology fixes the terminal set to {DAM, side_state}, so we set
    those terminal states MANUALLY on a GPCCA estimator over the PseudotimeKernel
    transition matrix and solve the SPARSE absorption-probability system (seconds
    on 26k cells) -- no dense Schur decomposition, no macrostate scan.

    The automatic macrostate-discovery path (compute_schur + compute_macrostates)
    is OFF by default (`auto=True` re-enables it): at 26k cells its real-Schur
    step densifies to an O(n^3) / 30 GB+ operation that dominated the entire
    build (~45 min) while only re-deriving terminal states the known biology
    already gives us. When enabled it now uses a SPARSE Krylov Schur and a single
    bounded n_states so it can never regress to the dense monster.

    Returns (dam_fate ndarray [NaN on failure], status string). Never raises."""
    n = adata.n_obs
    try:
        from cellrank.estimators import GPCCA
        from cellrank.kernels import PseudotimeKernel

        pk = PseudotimeKernel(adata, time_key="dpt_pseudotime")
        pk.compute_transition_matrix()
        g = GPCCA(pk)

        dam_fate = None
        status = None

        # --- optional automatic macrostates (sparse Krylov Schur, bounded) ---
        if auto:
            try:
                g.compute_schur(n_components=20, method="krylov")
                g.compute_macrostates(n_states=8, cluster_key="state")
                g.predict_terminal_states()
                g.compute_fate_probabilities()
                fp = g.fate_probabilities
                names = list(fp.names)
                dam_cols = [nm for nm in names
                            if nm.upper().startswith(dam_state.upper())]
                if dam_cols:
                    dam_fate = np.asarray(fp[:, dam_cols[0]].X).ravel()
                    status = f"auto macrostates; terminals={names}; dam_col={dam_cols[0]}"
            except Exception as e:  # noqa: BLE001
                status = f"auto failed ({type(e).__name__}: {str(e)[:120]})"

        # --- default: manual terminal set = {DAM, side_state} ----------------
        # Cheap sparse absorption solve. If the estimator demands a Schur basis
        # first, fall back to a SPARSE Krylov partial Schur (never dense Brandts).
        if dam_fate is None:
            term = pd.Series(pd.NA, index=adata.obs_names, dtype="object")
            st = adata.obs["state"].to_numpy()
            term[st == dam_state] = dam_state
            term[st == side_state] = side_state
            term = term.astype("category")
            try:
                g.set_terminal_states(term)
                g.compute_fate_probabilities()
            except Exception:  # noqa: BLE001 -- needs a Schur basis first
                g.compute_schur(n_components=10, method="krylov")
                g.set_terminal_states(term)
                g.compute_fate_probabilities()
            fp = g.fate_probabilities
            dam_fate = np.asarray(fp[:, dam_state].X).ravel()
            status = (status or "") + f" | manual terminal={{{dam_state},{side_state}}}"

        return dam_fate.astype(float), status
    except Exception as e:  # noqa: BLE001
        return np.full(n, np.nan), f"cellrank failed: {type(e).__name__}: {str(e)[:160]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--export-dir", default="storage/cache/trajectory")
    ap.add_argument("--n-neighbors", type=int, default=30)
    ap.add_argument("--n-dcs", type=int, default=10)
    ap.add_argument("--root-state", default="homeostatic")
    ap.add_argument("--dam-state", default="DAM")
    ap.add_argument("--side-state", default="IFN")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cellrank-auto", action="store_true",
                    help="re-enable CellRank automatic macrostate discovery "
                         "(sparse Krylov Schur); default uses the cheap manual "
                         "{DAM, side_state} terminal set")
    args = ap.parse_args()

    import scanpy as sc
    import anndata as ad

    sc.settings.verbosity = 1
    np.random.seed(args.seed)

    edir = args.export_dir
    emb = pd.read_csv(os.path.join(edir, "embedding.csv"), index_col=0)
    obs = pd.read_csv(os.path.join(edir, "obs.csv"), index_col=0)
    obs = obs.loc[emb.index]  # align order
    log(f"loaded embedding {emb.shape} + obs {obs.shape}")
    log(f"states: {obs['state'].value_counts().to_dict()}")

    X = emb.to_numpy(dtype=np.float32)
    adata = ad.AnnData(X=X.copy(), obs=obs.copy())
    adata.obsm["X_harmony"] = X
    adata.obs["state"] = adata.obs["state"].astype("category")

    # ---- neighbours + diffusion map on the SAME harmony space Slingshot uses
    sc.pp.neighbors(adata, use_rep="X_harmony", n_neighbors=args.n_neighbors,
                    random_state=args.seed)
    sc.tl.diffmap(adata, n_comps=max(15, args.n_dcs + 1))

    # ---- DPT from a principled homeostatic-tip root ----------------------
    iroot = pick_root_cell(adata, args.root_state, args.dam_state)
    adata.uns["iroot"] = iroot
    sc.tl.dpt(adata, n_dcs=args.n_dcs)
    log(f"DPT done; root cell idx={iroot} (state={adata.obs['state'].iloc[iroot]}); "
        f"pt range [{adata.obs['dpt_pseudotime'].min():.3f}, "
        f"{adata.obs['dpt_pseudotime'].max():.3f}]")

    # ---- checkpoint: DPT-only CSV so a later CellRank stall/crash can never
    # cost the (already-computed) DPT pseudotime. Overwritten with the full
    # table once CellRank returns.
    out_path = os.path.join(edir, "python_pseudotime.csv")
    pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "dpt_pseudotime": adata.obs["dpt_pseudotime"].to_numpy(),
            "cellrank_dam_fate": np.full(adata.n_obs, np.nan),
        }
    ).to_csv(out_path, index=False)
    log(f"checkpoint: wrote DPT-only {out_path}")

    # ---- PAGA topology (provenance: independent check of Slingshot's MST) -
    # Non-blocking: a PAGA failure must never cost the essential DPT/CellRank
    # per-cell outputs written below.
    state_order = list(adata.obs["state"].cat.categories)
    try:
        sc.tl.paga(adata, groups="state")
        paga_conn = np.asarray(adata.uns["paga"]["connectivities"].todense())
    except Exception as e:  # noqa: BLE001
        log(f"PAGA failed ({type(e).__name__}: {str(e)[:120]}); connectivities=null")
        paga_conn = None

    # ---- CellRank2 velocity-free DAM-fate probability (best-effort) -------
    dam_fate, cr_status = run_cellrank(adata, args.dam_state, args.side_state,
                                       args.seed, auto=args.cellrank_auto)
    log(f"CellRank: {cr_status}")
    n_finite = int(np.isfinite(dam_fate).sum())
    log(f"CellRank DAM-fate finite for {n_finite}/{adata.n_obs} cells")

    # ---- overwrite checkpoint with the full per-cell pseudotime CSV ------
    out = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "dpt_pseudotime": adata.obs["dpt_pseudotime"].to_numpy(),
            "cellrank_dam_fate": dam_fate,
        }
    )
    out.to_csv(out_path, index=False)
    log(f"wrote {out_path} ({out.shape[0]} cells)")

    # ---- provenance -------------------------------------------------------
    prov = {
        "n_cells": int(adata.n_obs),
        "n_dims_embedding": int(emb.shape[1]),
        "n_neighbors": args.n_neighbors,
        "n_dcs": args.n_dcs,
        "root_state": args.root_state,
        "root_cell_index": iroot,
        "root_cell_id": str(adata.obs_names[iroot]),
        "dpt_range": [float(adata.obs["dpt_pseudotime"].min()),
                      float(adata.obs["dpt_pseudotime"].max())],
        "paga_state_order": state_order,
        "paga_connectivities": (paga_conn.round(4).tolist()
                                if paga_conn is not None else None),
        "cellrank_status": cr_status,
        "cellrank_dam_fate_n_finite": n_finite,
        "versions": {
            "python": sys.version.split()[0],
            "scanpy": sc.__version__,
            "anndata": ad.__version__,
        },
    }
    try:
        import cellrank as cr
        prov["versions"]["cellrank"] = cr.__version__
    except Exception:  # noqa: BLE001
        prov["versions"]["cellrank"] = "unavailable"
    prov_path = os.path.join(edir, "python_provenance.json")
    with open(prov_path, "w") as fh:
        json.dump(prov, fh, indent=2)
    log(f"wrote {prov_path}")
    log("done.")


if __name__ == "__main__":
    main()
