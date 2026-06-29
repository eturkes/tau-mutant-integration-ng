#!/usr/bin/env python
"""scripts/build_scenic_grn.py

SCENIC arc (K2): the Python half of the bridge. Consumes the sparse bundle
written by scripts/export_microglia_for_scenic.R, builds a pyscenic-shaped
loom, and runs GRNBoost2 co-expression GRN inference under the locked
10-run consensus policy (anti-anchoring: run count + recurrence fixed in K1
BEFORE inspecting which TFs survive). See storage/notes/scenic_regulons_plan.md.

Why a separate Python script: pyscenic 0.12.1 needs Python <=3.10 / numpy
<1.24 and lives in the project-local micromamba env `scenic`, NOT the 3.12
.venv. R cannot drive it; this script runs IN the env.

GRN engine: the dask scheduler crashes on this node (TBB fork), so we use
arboreto's no-dask multiprocessing helper `arboreto_with_multiprocessing.py`.
LOAD-BEARING GOTCHA (K1): running that helper from its install dir
(.../site-packages/pyscenic/cli/) puts that dir on sys.path[0], where the
sibling module pyscenic.py shadows the pyscenic *package* ("'pyscenic' is not
a package"). FIX: copy the helper to a neutral dir (here, the scenic cache
dir) and invoke the copy. This script automates that copy.

Inputs (from export_microglia_for_scenic.R; storage/cache/scenic/):
  microglia_counts.mtx     genes x cells sparse integer counts (MatrixMarket)
  microglia_genes.txt      MGI symbols, row order
  microglia_cells.txt      cell IDs, col order
  microglia_colattrs.tsv   CellID, genotype, batch, genotype_batch
  ../../data/cistarget/allTFs_mm.txt   candidate regulators (1860 mouse TFs)

Outputs (storage/cache/scenic/):
  microglia.loom           genes x cells loom (Gene row attr; CellID + meta col attrs)
  adj_seed{S}.tsv          TF-target-importance adjacencies, one per seed (RESUMABLE)
  grn_progress.log         per-run wall time + running ETA

Run (in the env, backgrounded for the long multi-run):
  .micromamba/bin/micromamba run -n scenic \
      python scripts/build_scenic_grn.py --max-runs 10 --num-workers 8
The first invocation builds the loom and times run 1 (calibration); re-running
resumes (existing adj_seed*.tsv are skipped).
"""
import argparse
import os
import subprocess
import sys
import time

import loompy as lp
import numpy as np
import pandas as pd
import pyscenic
import scipy.io

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCENIC = os.path.join(ROOT, "storage", "cache", "scenic")
CISTARGET = os.path.join(ROOT, "storage", "data", "cistarget")

MTX = os.path.join(SCENIC, "microglia_counts.mtx")
GENES = os.path.join(SCENIC, "microglia_genes.txt")
CELLS = os.path.join(SCENIC, "microglia_cells.txt")
COLATTRS = os.path.join(SCENIC, "microglia_colattrs.tsv")
LOOM = os.path.join(SCENIC, "microglia.loom")
ALLTFS = os.path.join(CISTARGET, "allTFs_mm.txt")
PROGRESS = os.path.join(SCENIC, "grn_progress.log")

# Helper lives inside the installed pyscenic package; copy it out to dodge the
# sys.path[0] package-shadow gotcha (see module docstring).
HELPER_SRC = os.path.join(
    os.path.dirname(pyscenic.__file__), "cli", "arboreto_with_multiprocessing.py")
HELPER_COPY = os.path.join(SCENIC, "_grnboost2_mp.py")

# The stock helper dispatches one target at a time (chunksize=1), so with fast
# per-target regressions the 8 workers idle (~15% CPU) waiting on the parent's
# per-task round-trip -- dispatch-bound, ~5.5 cores wasted. Batching the imap
# dispatch is OUTPUT-INVARIANT (identical per-target regressions, same concat)
# and saturates the cores (~3x faster). See K2 PROGRESS NOTE.
CHUNKSIZE = 32


def prepare_helper():
    """Copy the helper to a neutral dir (dodge the package-shadow gotcha) and
    patch chunksize=1 -> CHUNKSIZE. Raise if the pattern is gone (don't silently
    run slow)."""
    src = open(HELPER_SRC).read()
    patched, n = src.replace("chunksize=1,", f"chunksize={CHUNKSIZE},"), \
        src.count("chunksize=1,")
    if n != 1:
        raise RuntimeError(
            f"expected exactly one 'chunksize=1,' in helper, found {n}; "
            "inspect arboreto_with_multiprocessing.py before patching.")
    with open(HELPER_COPY, "w") as fh:
        fh.write(patched)


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(PROGRESS, "a") as fh:
        fh.write(line + "\n")


def build_loom(rebuild=False):
    if os.path.exists(LOOM) and not rebuild:
        with lp.connect(LOOM, mode="r") as ds:
            log(f"loom exists: {LOOM} ({ds.shape[0]} genes x {ds.shape[1]} cells) -- skip build")
        return
    log("building loom from sparse bundle ...")
    mat = scipy.io.mmread(MTX).tocsc().astype("float32")   # genes x cells
    genes = np.array([ln.strip() for ln in open(GENES)])
    cells = np.array([ln.strip() for ln in open(CELLS)])
    attrs = pd.read_csv(COLATTRS, sep="\t", dtype=str)
    assert mat.shape == (len(genes), len(cells)), (mat.shape, len(genes), len(cells))
    assert list(attrs["CellID"]) == list(cells), "colattrs CellID order != cells.txt"
    row_attrs = {"Gene": genes}
    col_attrs = {
        "CellID": cells,
        "genotype": attrs["genotype"].to_numpy(),
        "batch": attrs["batch"].to_numpy(),
        "genotype_batch": attrs["genotype_batch"].to_numpy(),
    }
    lp.create(LOOM, mat, row_attrs, col_attrs)
    log(f"wrote loom: {LOOM} ({mat.shape[0]} genes x {mat.shape[1]} cells, nnz={mat.nnz})")


def run_grn(seeds, num_workers):
    prepare_helper()
    n_tfs = sum(1 for _ in open(ALLTFS))
    log(f"GRNBoost2: {len(seeds)} run(s) seeds={seeds}, num_workers={num_workers}, "
        f"{n_tfs} candidate TFs, chunksize={CHUNKSIZE}, helper={HELPER_COPY}")
    # Pin BLAS to 1 thread/worker so 8 processes do not oversubscribe 8 cores.
    env = dict(os.environ,
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    times = []
    for i, seed in enumerate(seeds, 1):
        out = os.path.join(SCENIC, f"adj_seed{seed}.tsv")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            n = sum(1 for _ in open(out)) - 1
            log(f"run {i}/{len(seeds)} seed={seed}: exists ({n} adjacencies) -- skip")
            continue
        log(f"run {i}/{len(seeds)} seed={seed}: starting -> {os.path.basename(out)}")
        t0 = time.time()
        # DENSE path (no --sparse): arboreto 0.1.6 core.py:125 calls `.A` on a
        # sparse target column, which scipy >=1.14 removed -> every target fails
        # under --sparse. Dense loads ~2.4 GB (26k x 11.5k float64), fine on 62 GB,
        # and is the path the K1 smoke validated.
        subprocess.run(
            [sys.executable, HELPER_COPY, LOOM, ALLTFS,
             "--method", "grnboost2", "--output", out,
             "--num_workers", str(num_workers), "--seed", str(seed)],
            env=env, check=True,
        )
        dt = time.time() - t0
        times.append(dt)
        n = sum(1 for _ in open(out)) - 1
        eta = (np.mean(times) * (len(seeds) - i)) if times else 0
        log(f"run {i}/{len(seeds)} seed={seed}: DONE in {dt/60:.1f} min, "
            f"{n} adjacencies; mean/run={np.mean(times)/60:.1f} min, "
            f"ETA remaining ~{eta/3600:.1f} h")
    log("GRNBoost2 orchestration complete for requested seeds.")


def main():
    ap = argparse.ArgumentParser(description="K2: build loom + run GRNBoost2 (10-run consensus).")
    ap.add_argument("--max-runs", type=int, default=10, help="number of seeds 1..N (default 10)")
    ap.add_argument("--num-workers", type=int, default=8, help="multiprocessing workers (default 8)")
    ap.add_argument("--rebuild-loom", action="store_true", help="force loom rebuild")
    ap.add_argument("--only-loom", action="store_true", help="build loom then exit")
    args = ap.parse_args()

    os.makedirs(SCENIC, exist_ok=True)
    build_loom(rebuild=args.rebuild_loom)
    if args.only_loom:
        log("--only-loom set; exiting after loom build.")
        return
    run_grn(seeds=list(range(1, args.max_runs + 1)), num_workers=args.num_workers)


if __name__ == "__main__":
    main()
