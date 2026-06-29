#!/usr/bin/env python
"""scripts/build_scenic_ctx_aucell.py  (SCENIC arc K3)

cisTarget motif pruning -> signed per-seed regulons -> 10-run >=80% (>=8/10)
recurrence consensus -> AUCell per-cell activity. Second compute half of the
SCENIC lane; consumes the GRNBoost2 adjacencies from K2 (build_scenic_grn.py).
See storage/notes/scenic_regulons_plan.md (step K3) for the locked policy.

WHY a self-gating background job. K2's 10-run GRNBoost2 (~7 h) may still be in
flight when this is launched. With --wait-for-adj this script POLLS until all 10
adj_seed*.tsv exist non-empty AND the GRN process (grn.pid) has exited, then runs
ctx with the full core budget (no contention). So the whole K2->K3 compute can be
kicked off once and left to finish unattended. RESUMABLE: existing per-seed
reg_seed*.csv are skipped, so a crash/relaunch costs only the unfinished seeds.

THE LOAD-BEARING DESIGN DECISION (consensus definition; fixed before inspecting
which TFs survive, per the K1 anti-anchoring lock "keep regulons recurrent in
>=80% of runs"):
  * Run `pyscenic ctx` independently on EACH of the 10 adjacency runs -> 10 signed
    regulon sets (canonical pySCENIC multi-run robustness protocol; ctx is the
    cheap, deterministic-given-its-input step).
  * EDGE-LEVEL recurrence is the consensus substrate (this is how SCENIC reports
    high-confidence links, Aibar 2017): a TF->target edge is high-confidence iff it
    recurs in >= RECURRENCE (=8) of the 10 runs. A consensus regulon = a TF(sign)
    with its high-confidence target set; kept iff it has >= MIN_TARGETS (=5,
    mirroring the section-14 decoupleR `min_targets=5L`) such edges. Because an
    edge can only reach 8/10 if its TF(sign) regulon existed in >=8 runs,
    regulon-level recurrence (>=80%) is implied -- edge-level is the stricter,
    target-set-quality-preserving form of the same locked rule.
  * RECOVERY CENSUS is reported separately at the looser REGULON-PRESENCE level
    (did TF X form a motif-supported regulon in >=8/10 runs at all?), because
    "is Spi1 recovered?" is a presence question, distinct from "what is its
    high-confidence target set?". Both are emitted, honestly, for every section-14
    / section-18 comparison TF.

Signed regulons: ctx is run with --all_modules + --expression_mtx_fname (the loom),
so targets are split by TF-target correlation into activating (+) and repressing
(-) regulons. This mirrors the SIGNED CollecTRI network used in section 14, making
the K4 head-to-head a clean "swap only the network" comparison. Caveat (documented
for K5): SCENIC repressing (-) regulons are less established than activating ones;
sign is reported throughout, never hidden. ctx motif-enrichment thresholds are the
pyscenic defaults (rank_threshold=5000, auc_threshold=0.05, nes_threshold=3.0,
min_genes=20) -- recorded in the progress log, applied as stated.

Outputs (storage/cache/scenic/):
  reg_seed{S}.csv                 per-seed ctx motif-enrichment table (RESUMABLE)
  scenic_consensus_regulons.gmt   consensus signed regulons (name=TF(+/-)) for AUCell
  scenic_regulons.tsv             long: TF,sign,target,recurrence,mean_importance (for decoupleR/K4)
  scenic_recurrence_dist.tsv      n edges / n regulons surviving at each threshold 1..10 (sensitivity diagnostic)
  scenic_recovery_census.tsv      section-14/18 comparison TFs: present? n runs? sign? n high-conf targets?
  aucell.csv.gz                   per-cell AUCell activity (cells x consensus regulons)
  ctx_progress.log                per-seed ctx wall time + consensus/aucell summary

Run (in the env, backgrounded; self-gates on K2):
  .micromamba/bin/micromamba run -n scenic \
      python scripts/build_scenic_ctx_aucell.py --wait-for-adj --num-workers 8
Validate the consensus logic with no env / no compute:
  python3 scripts/build_scenic_ctx_aucell.py --self-test
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCENIC = os.path.join(ROOT, "storage", "cache", "scenic")
CISTARGET = os.path.join(ROOT, "storage", "data", "cistarget")

LOOM = os.path.join(SCENIC, "microglia.loom")
RANKINGS = os.path.join(
    CISTARGET,
    "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIFS = os.path.join(CISTARGET, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
GRN_PID = os.path.join(SCENIC, "grn.pid")
PROGRESS = os.path.join(SCENIC, "ctx_progress.log")

GMT = os.path.join(SCENIC, "scenic_consensus_regulons.gmt")
LONG_TSV = os.path.join(SCENIC, "scenic_regulons.tsv")
DIST_TSV = os.path.join(SCENIC, "scenic_recurrence_dist.tsv")
CENSUS_TSV = os.path.join(SCENIC, "scenic_recovery_census.tsv")
AUCELL = os.path.join(SCENIC, "aucell.csv")
AUCELL_GZ = AUCELL + ".gz"

# section-14 verdict TFs + section-18 NF-kB family (the head-to-head targets,
# from the plan's "comparison ledger"). Sign-agnostic for the presence census.
CENSUS_TFS = {
    "amyloid_activation": ["Spi1", "Nfkb1", "Sp3"],
    "interaction_metabolic": ["Myc", "Creb1", "Tp53", "Jun"],
    "nfkb_family_sec18": ["Rela", "Nfkb1", "Nfkb2", "Relb", "Rel"],
}
NAME_RE = re.compile(r"^(.*?)\(([+-])\)$")  # "Spi1(+)" -> ("Spi1", "+")


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(PROGRESS, "a") as fh:
        fh.write(line + "\n")


# ----------------------------------------------------------------------------
# Consensus core (pure; no pyscenic dependency so --self-test runs under any
# python). Operates on per_seed: list (one per run) of {(tf, sign): {target: w}}.
# ----------------------------------------------------------------------------
def build_consensus(per_seed, recurrence, min_targets):
    """Edge-level >=`recurrence`/n-run consensus.

    Returns (consensus, edge_seedcount, reg_seedcount) where
      consensus[(tf, sign)] = {target: mean_importance over runs the edge appeared}
        for regulons with >= min_targets high-confidence edges, and
      edge_seedcount[(tf, sign, target)] = number of runs the edge appeared in,
      reg_seedcount[(tf, sign)]          = number of runs the regulon appeared in.
    """
    edge_seedcount = defaultdict(int)
    edge_weightsum = defaultdict(float)
    reg_seedcount = defaultdict(int)
    for seed_reg in per_seed:
        for (tf, sign), t2w in seed_reg.items():
            reg_seedcount[(tf, sign)] += 1
            for target, w in t2w.items():
                edge_seedcount[(tf, sign, target)] += 1
                edge_weightsum[(tf, sign, target)] += float(w)
    consensus = defaultdict(dict)
    for (tf, sign, target), c in edge_seedcount.items():
        if c >= recurrence:
            consensus[(tf, sign)][target] = edge_weightsum[(tf, sign, target)] / c
    consensus = {k: v for k, v in consensus.items() if len(v) >= min_targets}
    return consensus, dict(edge_seedcount), dict(reg_seedcount)


def recurrence_distribution(edge_seedcount, n_runs, min_targets):
    """For each threshold t in 1..n_runs: n high-confidence edges and n regulons
    that would retain >= min_targets edges. Sensitivity diagnostic; the headline
    uses the locked threshold only."""
    rows = []
    for t in range(1, n_runs + 1):
        n_edges = sum(1 for c in edge_seedcount.values() if c >= t)
        reg_tc = Counter((tf, sign)
                         for (tf, sign, _t), c in edge_seedcount.items() if c >= t)
        n_reg = sum(1 for v in reg_tc.values() if v >= min_targets)
        rows.append((t, n_edges, n_reg))
    return rows


def _self_test():
    # 3 runs, recurrence>=2, min_targets>=2.
    # A(+): t1 in all 3, t2 in 2, t3 in 1 -> kept {t1,t2}
    # B(-): t1,t2 in 2 each -> kept {t1,t2}
    # C(+): only t1 in 2 runs -> 1 edge < min_targets -> dropped (but present 2/3)
    per_seed = [
        {("A", "+"): {"t1": 5.0, "t2": 3.0, "t3": 1.0}, ("B", "-"): {"t1": 2.0, "t2": 2.0},
         ("C", "+"): {"t1": 1.0}},
        {("A", "+"): {"t1": 7.0, "t2": 3.0}, ("B", "-"): {"t1": 4.0, "t2": 2.0},
         ("C", "+"): {"t1": 1.0}},
        {("A", "+"): {"t1": 6.0}},
    ]
    cons, edge_sc, reg_sc = build_consensus(per_seed, recurrence=2, min_targets=2)
    assert set(cons.keys()) == {("A", "+"), ("B", "-")}, cons.keys()
    assert set(cons[("A", "+")]) == {"t1", "t2"}, cons[("A", "+")]
    assert abs(cons[("A", "+")]["t1"] - 18.0 / 3) < 1e-9, cons[("A", "+")]["t1"]
    assert set(cons[("B", "-")]) == {"t1", "t2"}
    assert reg_sc[("A", "+")] == 3 and reg_sc[("C", "+")] == 2
    assert edge_sc[("A", "+", "t1")] == 3 and edge_sc[("A", "+", "t3")] == 1
    dist = recurrence_distribution(edge_sc, n_runs=3, min_targets=2)
    assert dist[0][0] == 1 and dist[-1][0] == 3  # thresholds labelled 1..3
    # at t=3: only A(+)/t1 edge -> A has 1 edge < min_targets -> 0 regulons
    assert dist[2] == (3, 1, 0), dist[2]
    print("SELF-TEST PASSED: consensus, weighting, recurrence distribution.")


# ----------------------------------------------------------------------------
# pyscenic-backed steps (lazy imports so --self-test needs no env)
# ----------------------------------------------------------------------------
def _pyscenic_bin():
    return os.path.join(os.path.dirname(sys.executable), "pyscenic")


def _pid_alive(pidfile):
    if not os.path.exists(pidfile):
        return False
    try:
        pid = int(open(pidfile).read().strip())
    except (ValueError, OSError):
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def wait_for_adjacencies(seeds, poll, max_wait):
    """Block until every adj_seed{S}.tsv is non-empty, the GRN process has exited,
    and the youngest adjacency file is >=90 s old (fully flushed). Anti-contention
    gate so ctx never competes with GRNBoost2 for the 8 cores."""
    adj = [os.path.join(SCENIC, f"adj_seed{s}.tsv") for s in seeds]
    waited = 0
    while True:
        present = [p for p in adj if os.path.exists(p) and os.path.getsize(p) > 0]
        grn_running = _pid_alive(GRN_PID)
        young = max((time.time() - os.path.getmtime(p) for p in present), default=1e9)
        if len(present) == len(adj) and not grn_running and young >= 90:
            log(f"all {len(adj)} adjacency files present, GRN exited, files settled "
                f"(youngest {young:.0f}s old) -- proceeding to ctx.")
            return
        log(f"waiting for K2: {len(present)}/{len(adj)} adj files, "
            f"GRN_running={grn_running}, youngest={young:.0f}s; sleep {poll}s "
            f"(waited {waited/60:.0f} min).")
        if waited >= max_wait:
            raise TimeoutError(
                f"adjacencies not ready after {max_wait/3600:.1f} h "
                f"({len(present)}/{len(adj)} present, GRN_running={grn_running}).")
        time.sleep(poll)
        waited += poll


def run_ctx_per_seed(seeds, num_workers, all_modules, mask_dropouts):
    """`pyscenic ctx` per adjacency run -> reg_seed{S}.csv (RESUMABLE)."""
    pys = _pyscenic_bin()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    log(f"ctx: {len(seeds)} run(s), workers={num_workers}, all_modules={all_modules}, "
        f"mask_dropouts={mask_dropouts}, defaults rank=5000/auc=0.05/nes=3.0/min_genes=20.")
    times = []
    for i, seed in enumerate(seeds, 1):
        adj = os.path.join(SCENIC, f"adj_seed{seed}.tsv")
        out = os.path.join(SCENIC, f"reg_seed{seed}.csv")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            log(f"ctx {i}/{len(seeds)} seed={seed}: reg exists -- skip")
            continue
        cmd = [pys, "ctx", adj, RANKINGS,
               "--annotations_fname", MOTIFS,
               "--expression_mtx_fname", LOOM,
               "--mode", "custom_multiprocessing",
               "--num_workers", str(num_workers),
               "--output", out]
        if all_modules:
            cmd.append("--all_modules")
        if mask_dropouts:
            cmd.append("--mask_dropouts")
        log(f"ctx {i}/{len(seeds)} seed={seed}: starting -> {os.path.basename(out)}")
        t0 = time.time()
        subprocess.run(cmd, env=env, check=True)
        dt = time.time() - t0
        times.append(dt)
        eta = (sum(times) / len(times)) * (len(seeds) - i)
        log(f"ctx {i}/{len(seeds)} seed={seed}: DONE in {dt/60:.1f} min; "
            f"mean/run={sum(times)/len(times)/60:.1f} min, ETA ~{eta/60:.0f} min")


def load_seed_regulons(seed):
    """reg_seed{S}.csv -> {(tf, sign): {target: weight}} via pyscenic df2regulons."""
    from pyscenic.transform import df2regulons
    from pyscenic.utils import load_motifs
    out = os.path.join(SCENIC, f"reg_seed{seed}.csv")
    regs = df2regulons(load_motifs(out))
    d = {}
    for r in regs:
        m = NAME_RE.match(r.name)
        if m:
            tf, sign = m.group(1), m.group(2)
        else:  # unsigned -> treat as activating, note it
            tf, sign = r.transcription_factor, "+"
        d[(tf, sign)] = dict(r.gene2weight)
    return d


def write_gmt(consensus, edge_seedcount, path):
    with open(path, "w") as fh:
        for (tf, sign), t2w in sorted(consensus.items()):
            name = f"{tf}({sign})"
            rec = min(edge_seedcount[(tf, sign, t)] for t in t2w)  # min recurrence in set
            desc = f"scenic_consensus_minrec{rec}_n{len(t2w)}"
            fh.write("\t".join([name, desc, *sorted(t2w)]) + "\n")


def write_long_tsv(consensus, edge_seedcount, path):
    with open(path, "w") as fh:
        fh.write("TF\tsign\ttarget\trecurrence\tmean_importance\n")
        for (tf, sign), t2w in sorted(consensus.items()):
            for target in sorted(t2w):
                fh.write(f"{tf}\t{sign}\t{target}\t"
                         f"{edge_seedcount[(tf, sign, target)]}\t{t2w[target]:.6g}\n")


def write_recurrence_dist(edge_seedcount, n_runs, min_targets, path):
    rows = recurrence_distribution(edge_seedcount, n_runs, min_targets)
    with open(path, "w") as fh:
        fh.write("min_runs\tn_highconf_edges\tn_regulons_ge_min_targets\n")
        for t, ne, nr in rows:
            fh.write(f"{t}\t{ne}\t{nr}\n")
    return rows


def write_recovery_census(consensus, reg_seedcount, recurrence, path):
    """Per comparison TF: best regulon-presence (n runs) over signs, whether it
    clears the >=recurrence presence bar, and its consensus high-confidence target
    count. Reported sign-agnostic for presence, sign-resolved for the kept regulon."""
    cons_by_tf = defaultdict(list)  # tf -> [(sign, n_targets)]
    for (tf, sign), t2w in consensus.items():
        cons_by_tf[tf].append((sign, len(t2w)))
    pres_by_tf = defaultdict(list)  # tf -> [(sign, n_runs)]
    for (tf, sign), n in reg_seedcount.items():
        pres_by_tf[tf].append((sign, n))
    with open(path, "w") as fh:
        fh.write("axis\tTF\tregulon_present_ge_thr\tmax_runs_present\tpresent_signs"
                 "\tin_consensus\tconsensus_signs_ntargets\n")
        for axis, tfs in CENSUS_TFS.items():
            for tf in tfs:
                pres = pres_by_tf.get(tf, [])
                max_runs = max((n for _s, n in pres), default=0)
                present_signs = ",".join(f"{s}:{n}" for s, n in sorted(pres)) or "-"
                in_cons = tf in cons_by_tf
                cons_str = ",".join(f"{s}:{n}" for s, n in sorted(cons_by_tf.get(tf, []))) or "-"
                fh.write(f"{axis}\t{tf}\t{int(max_runs >= recurrence)}\t{max_runs}"
                         f"\t{present_signs}\t{int(in_cons)}\t{cons_str}\n")


def run_aucell(num_workers):
    """AUCell over the consensus regulon gmt -> cells x regulons, gzipped."""
    pys = _pyscenic_bin()
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    log(f"aucell: scoring {LOOM} against {os.path.basename(GMT)} ...")
    t0 = time.time()
    subprocess.run([pys, "aucell", LOOM, GMT, "--output", AUCELL,
                    "--num_workers", str(num_workers)], env=env, check=True)
    with open(AUCELL, "rb") as fi, gzip.open(AUCELL_GZ, "wb") as fo:
        fo.writelines(fi)
    os.remove(AUCELL)
    log(f"aucell DONE in {(time.time()-t0)/60:.1f} min -> {os.path.basename(AUCELL_GZ)}")


def main():
    ap = argparse.ArgumentParser(description="K3: ctx -> consensus regulons -> AUCell.")
    ap.add_argument("--max-runs", type=int, default=10, help="seeds 1..N (default 10)")
    ap.add_argument("--num-workers", type=int, default=8)
    ap.add_argument("--recurrence", type=int, default=8,
                    help="min runs an edge must appear in (locked >=8/10)")
    ap.add_argument("--min-targets", type=int, default=5,
                    help="min high-confidence targets to keep a regulon (mirror section 14)")
    ap.add_argument("--wait-for-adj", action="store_true",
                    help="poll until K2's adjacencies are ready + GRN exited")
    ap.add_argument("--poll", type=int, default=300, help="wait poll seconds (default 300)")
    ap.add_argument("--max-wait", type=float, default=12 * 3600, help="give up after (s)")
    ap.add_argument("--activating-only", action="store_true",
                    help="ctx without --all_modules (drop repressing regulons)")
    ap.add_argument("--mask-dropouts", action="store_true",
                    help="ctx --mask_dropouts (default off: zeros kept in correlation)")
    ap.add_argument("--skip-aucell", action="store_true", help="stop after consensus export")
    ap.add_argument("--self-test", action="store_true",
                    help="run the consensus unit test (no env / no compute) and exit")
    args = ap.parse_args()

    if args.self_test:
        _self_test()
        return

    os.makedirs(SCENIC, exist_ok=True)
    seeds = list(range(1, args.max_runs + 1))
    for p in (LOOM, RANKINGS, MOTIFS):
        if not os.path.exists(p):
            raise FileNotFoundError(p)

    if args.wait_for_adj:
        wait_for_adjacencies(seeds, poll=args.poll, max_wait=args.max_wait)

    run_ctx_per_seed(seeds, args.num_workers,
                     all_modules=not args.activating_only,
                     mask_dropouts=args.mask_dropouts)

    log("loading per-seed regulons + building consensus ...")
    per_seed = [load_seed_regulons(s) for s in seeds]
    for s, d in zip(seeds, per_seed):
        log(f"  seed {s}: {len(d)} signed regulons")
    consensus, edge_sc, reg_sc = build_consensus(
        per_seed, recurrence=args.recurrence, min_targets=args.min_targets)
    n_pos = sum(1 for (_t, sgn) in consensus if sgn == "+")
    n_neg = len(consensus) - n_pos
    log(f"CONSENSUS (>= {args.recurrence}/{len(seeds)} runs, >= {args.min_targets} "
        f"targets): {len(consensus)} regulons ({n_pos} activating, {n_neg} repressing).")
    if len(consensus) < 20:
        log(f"WARNING: only {len(consensus)} consensus regulons (<20). Low recovery "
            "is a known SCENIC single-cell-type failure mode -- inspect the "
            "recurrence distribution before interpreting (K3 verify note).")

    write_gmt(consensus, edge_sc, GMT)
    write_long_tsv(consensus, edge_sc, LONG_TSV)
    dist = write_recurrence_dist(edge_sc, len(seeds), args.min_targets, DIST_TSV)
    write_recovery_census(consensus, reg_sc, args.recurrence, CENSUS_TSV)
    log("recurrence sensitivity (min_runs: n_edges, n_regulons): "
        + "; ".join(f"{t}:{ne},{nr}" for t, ne, nr in dist))
    log(f"wrote {os.path.basename(GMT)}, {os.path.basename(LONG_TSV)}, "
        f"{os.path.basename(DIST_TSV)}, {os.path.basename(CENSUS_TSV)}")

    if args.skip_aucell:
        log("--skip-aucell set; stopping after consensus export.")
        return
    run_aucell(args.num_workers)
    log("K3 compute complete: consensus regulons + AUCell written. Mark K3 DONE "
        "after eyeballing scenic_recovery_census.tsv + scenic_recurrence_dist.tsv.")


if __name__ == "__main__":
    main()
