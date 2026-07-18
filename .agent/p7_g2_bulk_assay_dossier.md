# P7.2 G2 dossier — 24M bulk assay identity and run-order contract

## Decision

The historical `proteome` layer is **not a global-proteome measurement**. It is a
TiO2-enriched phospho-PTM Spectronaut report whose positive PTM-row intensities are
summed to protein groups. The historical `proteome_*` / `Proteome` code tokens and
list keys remain stable for pipeline compatibility; they are not an assay-identity
claim.

Both retained bulk layers come from the same `Naoto-Hippo_TiO2_DIA` phospho
acquisition:

- the historical `proteome` layer is the protein-group-sum view; and
- the `phospho` layer is the phosphosite-row view.

No separate global-proteome input exists under `storage/data/`.

## File and column evidence

- `storage/data/proteomics_nonfiltered_nonnormalised.tsv` has 30 columns, including
  phospho-PTM annotations `Gene-pSite`, `PTM.ProteinId`, `PTM.CollapseKey`,
  `PTM.Multiplicity`, `PTM.ModificationTitle`, `PTM.SiteAA`, `PTM.SiteLocation`,
  `PTM.FlankingRegion`, and `Phosphopeptide`, plus 16
  `Naoto-Hippo_TiO2_DIA_NN.raw.PTM.Quantity` intensity columns.
- `storage/data/phosphoproteomics_nonfiltered_nonnormalised.tsv` has 81 columns and
  carries the same `Naoto-Hippo_TiO2_DIA` acquisition at phosphosite level: 67
  `.PTM.Quantity` manifest columns exist, of which the 16-run 24M subset is used.
- TiO2 denotes titanium-dioxide phosphopeptide enrichment. These two files are two
  reporting levels of that phospho assay, not independent global-proteome and
  phosphoproteome assays.

## Code-path evidence

- `R/analysis/modality_de.R::match_24m_bulk_columns()` matches the exact 16-run 24M
  sample-key subset, orders it by the key, and assigns `run_index = 1:16`.
- `protein_group_features()` uses `PG.ProteinGroups` as the protein-group feature ID.
- `aggregate_proteome_raw()` converts non-positive values to missing, then sums the
  remaining positive `.raw.PTM.Quantity` PTM-row intensities by
  `PG.ProteinGroups` before log2 transformation.
- `prepare_proteome_24m_matrix()` performs the protein-group aggregation, log2
  transform, sample-wise median shift, and prevalence filter.
- `run_proteome_de_24m()` retains the historical name but records
  `feature_id = "PG.ProteinGroups"` and aggregation of positive peptide/PTM-row
  intensities; its figure-facing `$top` tables remain the primary no-batch fit.
- `run_phospho_de_24m()` applies the same 16-run match to the phosphosite-level
  report and fits the same five canonical contrasts.

**STANDING external-record request (user-supplied; not in-repo):** vendor TiO2
enrichment/fractionation records plus the acquisition method/template for the
`Naoto-Hippo_TiO2_DIA` run set. This does not block P7.2.

## Run-order contract and confound

The 24M acquisition order is perfectly genotype-blocked:

| Runs | Genotype |
|---|---|
| 01–04 | MAPTKI |
| 05–08 | P301S |
| 09–12 | NLGF_MAPTKI |
| 13–16 | NLGF_P301S |

The primary bulk design is `factorial_design(meta, add_batch = FALSE)`. Run order is
therefore a deterministic function of genotype at the between-block level, so the
between-genotype contrasts cannot be separated from acquisition order/batch.

P7.2 adds an integrity-record sensitivity only: the primary four-column factorial
design is augmented with one continuous mean-centred `run_index` column
(`-7.5, -6.5, ..., 7.5`). The augmented design has rank 5 and 11 residual degrees of
freedom. The same five contrast vectors are refit with zero weight on `run_index`.
Figures continue to read the unchanged primary `$top` logFCs.

The summaries below are signed feature-mean logFCs over features finite in both fits;
`max |feature shift|` is the largest absolute feature-level logFC change. The
continuous term captures only within-acquisition linear drift and does **not** repair
the genotype/order alias.

### Historical `proteome` layer — TiO2 phospho protein-group sums

| Contrast | Primary mean logFC | Run-index-adjusted mean logFC | Mean shift | Max \|feature shift\| |
|---|---:|---:|---:|---:|
| `tau_alone` | -0.081319760 | 0.040351662 | +0.121671422 | 8.139318051 |
| `nlgf_in_maptki` | 0.112308386 | 0.354895758 | +0.242587371 | 14.650772491 |
| `nlgf_in_p301s` | 0.024661037 | 0.268487994 | +0.243826956 | 14.325199769 |
| `tau_in_nlgf` | -0.168967109 | -0.046056101 | +0.122911008 | 7.813745329 |
| `interaction` | -0.087647349 | -0.086407764 | +0.001239585 | 2.027386897 |

### Phosphosite layer

| Contrast | Primary mean logFC | Run-index-adjusted mean logFC | Mean shift | Max \|feature shift\| |
|---|---:|---:|---:|---:|
| `tau_alone` | -0.141921341 | -0.184412971 | -0.042491630 | 11.388507948 |
| `nlgf_in_maptki` | -0.008073148 | -0.093595010 | -0.085521862 | 22.517637180 |
| `nlgf_in_p301s` | 0.087189662 | 0.002995360 | -0.084194302 | 23.232556213 |
| `tau_in_nlgf` | -0.046658530 | -0.087822601 | -0.041164071 | 12.299588584 |
| `interaction` | 0.095262811 | 0.096590370 | +0.001327559 | 4.426943032 |

The aggregate interaction shifts are approximately zero, as expected because the
2x2 interaction contrast is orthogonal to a linear run index for complete 16-run
profiles. Feature-specific missingness breaks that exact per-feature geometry:
2,478/3,379 protein-group rows and 9,261/17,707 phosphosite rows are complete, and
the maximum interaction shift among those complete rows is zero; the non-zero maxima
above arise from partially observed features. This sensitivity remains descriptive
and context-only.
