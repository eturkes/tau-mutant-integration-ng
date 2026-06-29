#!/usr/bin/env python3
"""scripts/build_human_microglia.py

Session H1 of the human cross-species validation plan
(storage/notes/human_validation_plan.md).

Acquires the SEA-AD (Gabitto 2024) microglia/PVM single-nucleus data for
the two cortical regions chosen at the H1 gate -- MTG and DLPFC -- and
emits the durable H1 artifacts: a per-donor neuropathology table with
NUMERIC amyloid (Thal) and tau (Braak) axes, plus a provenance record.

The source h5ads are the CELLxGENE Discover per-cell-type-per-region
splits (microglia-only, ~0.5-0.7 GB each), so NO subsetting of the 50 GB
"Whole Taxonomy" files is needed. They are downloaded once (by the H1
session, via curl) to storage/data/seaad/ and read here in backed mode
(only obs is materialised; the expression matrix stays on disk for H3
scoring).

Why ordinal stages, not continuous IHC: the CELLxGENE obs carries Thal
phase, Braak stage, CERAD, ADNC and the Continuous Pseudo-progression
Score, but NOT the per-donor image-quantified Abeta/AT8 densities (those
live in the SEA-AD-native supplements at brain-map.org). The amyloid axis
is therefore Thal phase (0-5) and the tau axis Braak stage (0-VI); CERAD
is a secondary amyloid (neuritic-plaque) axis and CPS a collapsed
severity axis. Joining the continuous IHC densities is a documented H1
enhancement, not a blocker.

Inputs (downloaded by the H1 session):
  storage/data/seaad/microglia_mtg.h5ad     (CELLxGENE dataset c66e3198-...)
  storage/data/seaad/microglia_dlpfc.h5ad    (CELLxGENE dataset a3198428-...)

Outputs:
  storage/cache/human_seaad_donor_neuropath.csv   per-donor x region table
  storage/cache/human_microglia_provenance.txt     full provenance + checks

Run: ./.venv/bin/python scripts/build_human_microglia.py
"""
import os
import sys
import textwrap
import numpy as np
import pandas as pd
import anndata as ad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEAAD = os.path.join(ROOT, "storage", "data", "seaad")
CACHE = os.path.join(ROOT, "storage", "cache")

# SEA-AD CELLxGENE Discover provenance (collection 1ca90a2d-...).
SOURCES = {
    "MTG": dict(
        path=os.path.join(SEAAD, "microglia_mtg.h5ad"),
        dataset_id="c76098ba-eed3-45b1-98f2-96fcac55ed18",
        url="https://datasets.cellxgene.cziscience.com/c66e3198-c766-499e-a609-7b462a41295b.h5ad",
    ),
    "DLPFC": dict(
        path=os.path.join(SEAAD, "microglia_dlpfc.h5ad"),
        dataset_id="100c6145-7b0e-4ba6-81c1-ffebed0d1ac4",
        url="https://datasets.cellxgene.cziscience.com/a3198428-dd16-4344-8564-3897c1ccdb4b.h5ad",
    ),
}

# --- ordinal neuropathology -> numeric axis recodings -----------------------
# "Reference" = Allen neurotypical reference donors (separate from the SEA-AD
# AD-continuum cohort); they carry no staged neuropath, so they map to NaN and
# are flagged is_reference for downstream filtering at H4.
THAL = {f"Thal {i}": i for i in range(6)}                       # amyloid 0-5
BRAAK = {"Braak 0": 0, "Braak II": 2, "Braak III": 3,           # tau (no Braak I)
         "Braak IV": 4, "Braak V": 5, "Braak VI": 6}
CERAD = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}  # neuritic plaque
ADNC = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}
COG = {"No dementia": 0, "Dementia": 1}
APOE4 = {"N": 0, "Y": 1}

DONOR_COLS = ["donor_id", "Thal phase", "Braak stage", "CERAD score", "ADNC",
              "Continuous Pseudo-progression Score", "Cognitive status",
              "APOE4 status", "Age at death", "PMI", "sex", "disease"]


def num(series, mapping):
    return series.astype(str).map(mapping).astype("float")


def age_num(series):
    # SEA-AD censors high ages as "90+ years"; pull the leading integer.
    return series.astype(str).str.extract(r"(\d+)").iloc[:, 0].astype("float")


def main():
    prov = []
    donor_frames = []
    var_ref = None
    supertype_counts = {}
    xtabs = {}

    for region, meta in SOURCES.items():
        if not os.path.exists(meta["path"]):
            sys.exit(f"MISSING input h5ad: {meta['path']} (download it first)")
        a = ad.read_h5ad(meta["path"], backed="r")
        o = a.obs
        prov.append(f"[{region}] {a.shape[0]} cells x {a.shape[1]} genes | "
                    f"{o['donor_id'].nunique()} donors | dataset {meta['dataset_id']}")

        # gene-space consistency across regions (H3 scores uniformly)
        vn = list(a.var["feature_name"].astype(str))
        if var_ref is None:
            var_ref = vn
        else:
            same = (len(vn) == len(var_ref)) and (vn == var_ref)
            prov.append(f"[{region}] gene space identical to first region: {same}")

        supertype_counts[region] = o["Supertype"].value_counts().to_dict()

        # Per-region obs differ slightly (DLPFC lacks the MTG-only Continuous
        # Pseudo-progression Score); select what exists, NaN-fill the rest.
        present = [c for c in DONOR_COLS if c in o.columns]
        missing = [c for c in DONOR_COLS if c not in o.columns]
        if missing:
            prov.append(f"[{region}] obs missing {missing} -> NaN-filled")
        d = o.drop_duplicates("donor_id")[present].copy()
        for c in missing:
            d[c] = np.nan
        d["region"] = region
        d["amyloid_thal"] = num(d["Thal phase"], THAL)
        d["tau_braak"] = num(d["Braak stage"], BRAAK)
        d["cerad"] = num(d["CERAD score"], CERAD)
        d["adnc"] = num(d["ADNC"], ADNC)
        d["cps"] = pd.to_numeric(d["Continuous Pseudo-progression Score"], errors="coerce")
        d["cognitive"] = num(d["Cognitive status"], COG)
        d["apoe4"] = num(d["APOE4 status"], APOE4)
        d["age"] = age_num(d["Age at death"])
        d["pmi"] = pd.to_numeric(d["PMI"].astype(str), errors="coerce")
        d["is_reference"] = d["Thal phase"].astype(str).eq("Reference")
        donor_frames.append(d)

        # donor-level identifiability cross-tab + collinearity (continuum only)
        cont = d[~d["is_reference"]]
        xt = pd.crosstab(cont["amyloid_thal"], cont["tau_braak"])
        xtabs[region] = xt
        rho = cont[["amyloid_thal", "tau_braak"]].corr(method="spearman").iloc[0, 1]
        prov.append(f"[{region}] continuum donors={len(cont)} | "
                    f"Spearman(amyloid_thal, tau_braak)={rho:.3f} "
                    f"(collinearity baseline; off-diagonal donors are the "
                    f"identifiability budget for the interaction term)")

    donors = pd.concat(donor_frames, ignore_index=True)
    out_csv = os.path.join(CACHE, "human_seaad_donor_neuropath.csv")
    donors.to_csv(out_csv, index=False)

    # ---- provenance file ----
    lines = []
    lines.append("SEA-AD microglia/PVM acquisition -- human cross-species validation H1")
    lines.append("=" * 72)
    lines.append("Source: CELLxGENE Discover collection "
                 "1ca90a2d-2943-483d-b678-b809bf464c30 (SEA-AD, Gabitto 2024).")
    for region, meta in SOURCES.items():
        lines.append(f"  {region}: {meta['url']}")
    lines.append("")
    lines.append("Per-region summary:")
    lines += [f"  {p}" for p in prov]
    lines.append("")
    lines.append("Amyloid axis = Thal phase (0-5); tau axis = Braak stage "
                 "(0,2-6; no Braak I in SEA-AD). CERAD (0-3) = neuritic-plaque "
                 "amyloid; ADNC (0-3) = combined ABC; CPS = continuous severity.")
    lines.append("Reference (Allen neurotypical) donors -> NaN axes, "
                 "is_reference=True; excluded from interaction modelling at H4.")
    lines.append("Expression: X = normalised log1p (CELLxGENE); raw counts in "
                 ".raw; gene symbols in var.feature_name. h5ads ARE the compact "
                 "microglia-only cache (no re-subset needed).")
    lines.append("")
    for region, xt in xtabs.items():
        lines.append(f"Donor-level amyloid_thal x tau_braak ({region}, continuum):")
        lines.append(textwrap.indent(xt.to_string(), "  "))
        lines.append("")
    for region, sc in supertype_counts.items():
        lines.append(f"Supertype (microglial state) counts ({region}):")
        for k, v in sc.items():
            lines.append(f"  {k}: {v}")
        lines.append("")
    prov_path = os.path.join(CACHE, "human_microglia_provenance.txt")
    with open(prov_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # ---- smoke summary ----
    print("WROTE", out_csv, f"({len(donors)} donor x region rows)")
    print("WROTE", prov_path)
    print("\nper-donor table columns:", list(donors.columns))
    print("\ncontinuum donors per region:")
    print(donors[~donors.is_reference].groupby("region")["donor_id"].nunique())
    print("\nneuropath non-null (continuum, MTG):")
    m = donors[(donors.region == "MTG") & (~donors.is_reference)]
    print(m[["amyloid_thal", "tau_braak", "cerad", "adnc", "cps", "cognitive",
             "apoe4", "age"]].notna().sum())


if __name__ == "__main__":
    main()
