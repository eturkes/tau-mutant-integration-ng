# P7.3 G3 dossier — source provenance and cross-modality claim gate

## Decision

Live report inputs = 16/67 bulk runs + 91/112 ideal-grid GeoMx AOIs + 16/16
snRNAseq replicate units. In-repo provenance resolves each modality's retained design but
contains **no cross-assay animal/aliquot crosswalk**. Animal-level or paired
cross-modality integration is therefore **PROHIBITED** until the user supplies and
validates that crosswalk.

Current report = compliant: cross-assay harmonization occurs only by gene symbol in the
four-method amyloid-response logFC scatter (Figure 8); retained bulk/GeoMx figures are
modality-native context. No within-animal or paired cross-modality claim is made.

## 24M bulk — 67-run manifest

`storage/data/proteomics_sample_key.csv` = 2 columns (`File name`,
`Sample/Condtion`; source misspelling retained) × 67 data rows. File names encode two
acquisition sets with pattern `..._Hippo_TiO2_DIA_NN.PTM.Quantity` (`NN` =
zero-padded run index): Naoto = 26 runs, rows 1–26; Set6 = 41 runs, rows 27–67.

### Naoto set — 26 runs

| Manifest rows | File run | Raw label | Count | Report disposition |
|---|---|---|---:|---|
| 1–4 | `Naoto-Hippo_TiO2_DIA_01..04` | `MAPT-KI_24M` | 4 | USED |
| 5–8 | `Naoto-Hippo_TiO2_DIA_05..08` | `P301S+3_24M` | 4 | USED |
| 9–12 | `Naoto-Hippo_TiO2_DIA_09..12` | `NLGF-MAPT-KI_24M` | 4 | USED |
| 13–16 | `Naoto-Hippo_TiO2_DIA_13..16` | `NLGF-P301S+3_24M` | 4 | USED |
| 17–20 | `Naoto-Hippo_TiO2_DIA_17..20` | `MAPT-KI_30M` | 4 | UNUSED |
| 21–26 | `Naoto-Hippo_TiO2_DIA_21..26` | `P301S+3_30M` | 6 | UNUSED |
| **Subtotal** |  |  | **26** | **16 used; 10 unused** |

Rows 1–16 = complete balanced 24M 2×2, four runs/genotype, and the **only runs
used in the report**. Rows 17–26 = 30M, unused and incomplete/unbalanced: 4
`MAPT-KI_30M` versus 6 `P301S+3_30M`, with no NLGF 30M arms. The 30M subset is
not a 2×2. The 24M block being the Naoto set's only complete balanced 2×2 is a
plausible selection rationale; the manifest does not itself record that rationale.

### Set6 set — 41 runs, all unused

| Manifest rows | Raw label | Count |
|---|---|---:|
| 27–29 | `MAPT-wt-H1-haplotype` | 3 |
| 30–32 | `MAPT-10+16` | 3 |
| 33–35 | `MAPT-N279K` | 3 |
| 36–40 | `MAPT-P301S-10+3` | 5 |
| 41–44 | `MAPT-S305N-10+3` | 4 |
| 45–48 | `MAPT-WT` | 4 |
| 49–51 | `MAPT-10+3` | 3 |
| 52–55 | `NLGF-MAPT-KI` | 4 |
| 56–59 | `NLGF-10+3` | 4 |
| 60–63 | `NLGF-S305N-10+3` | 4 |
| 64–67 | `NLGF-P301S-10+3` | 4 |
| **Subtotal** |  | **41** |

Manifest arithmetic: 26 Naoto + 41 Set6 = 67 runs.

### Set6 series attribution — literature-inferred, not manifest-certain

This attribution follows the P7.1 construct literature work; it is an inference from
label nomenclature and published model series, not animal-level provenance encoded by
the manifest.

- Benzow/Koob H1 `MAPT-GR` series — distinct human-tau targeted gene-replacement
  line (human MAPT genomic region replaces the mouse locus, not a random-insertion
  transgenic), H1 haplotype:
  `MAPT-wt-H1-haplotype` (matched WT control), `MAPT-10+16`, `MAPT-N279K`.
- Saito humanized base-edited MAPT-KI series — same platform as the report's
  MAPTKI/P301S groups: `MAPT-10+3`, `MAPT-P301S-10+3`,
  `MAPT-S305N-10+3`; `NLGF-MAPT-KI`, `NLGF-10+3`,
  `NLGF-S305N-10+3`, and `NLGF-P301S-10+3` denote App^NL-G-F amyloid
  crosses of that series.
- `MAPT-WT` (rows 45–48) = ambiguous from the manifest alone: potentially the Saito
  humanized WT-KI or a plain WT. User records must disambiguate it and confirm every
  series assignment at animal level.

**HARD CAUTION:** Set6 WT-like labels `MAPT-WT` and
`MAPT-wt-H1-haplotype` are not interchangeable with the Naoto report
`MAPT-KI`. They must not be pooled as a shared control.

## 16-of-67 bulk selection mechanism

`R/core/io.R::proteomics_sample_meta(n_keep = 16L)` makes the report selection
positionally, then validates its identity:

1. `readr::read_csv()` reads the 67-row manifest.
2. `dplyr::slice_head(n = 16)` takes rows 1–16 = the Naoto 24M block.
3. Four raw labels map to canonical genotypes `{MAPTKI, P301S, NLGF_MAPTKI,
   NLGF_P301S}`.
4. Hard `stopifnot` contract requires `nrow == 16`, no missing genotype,
   `setequal(label, {MAPT-KI_24M, P301S+3_24M, NLGF-MAPT-KI_24M,
   NLGF-P301S+3_24M})`, exactly 4 runs/genotype, and unique `file_name` +
   `col_stub`.

Selection is positional but fail-loud validated to be exactly the balanced 24M 2×2.
Naoto 30M rows 17–26 and all Set6 rows 27–67 never reach any target.

## GeoMx — 91-of-112 AOI gap

`storage/data/geomx.rds` = processed Seurat object with 91 AOIs. Design metadata
includes `genotype`, `bio_rep`, `tech_rep`, `slide_rep`, `roi`, `segment`, `aoi`,
`area`, `nuclei`, `SampleID`, plus `slide name` / `Scan Name` physical-DSP-scan
fields and WTA negative-probe/q-normalization QC.

- `bio_rep == slide_rep`; distributions are identical.
- `tech_rep` levels = 1–7.
- One segment only: `segment = "Segment 1"`; `aoi = "Segment 1-aoi-001"` for all
  91 AOIs.
- `SampleID` = 91 unique DSP barcodes across at least two DSP scan plates:
  `DSP-1001660019825-A` and `DSP-1001660022195-B`.

Ideal grid = 4 genotypes × 4 `bio_rep` × 7 `tech_rep` = 112 AOIs. Observed = 91;
**21 AOIs are missing across 7 non-complete (`genotype`, `bio_rep`) cells**. The
other 9 cells are complete at 7 AOIs.

| Genotype | `bio_rep` | Observed AOIs | Present `tech_rep` | Missing `tech_rep` | Deficit |
|---|---:|---:|---|---|---:|
| `MAPTKI` | 1 | 0 | `{}` | `{1,2,3,4,5,6,7}` | −7 |
| `MAPTKI` | 4 | 6 | `{1,2,4,5,6,7}` | `{3}` | −1 |
| `P301S` | 1 | 2 | `{6,7}` | `{1,2,3,4,5}` | −5 |
| `P301S` | 2 | 5 | `{1,2,3,5,6}` | `{4,7}` | −2 |
| `P301S` | 4 | 6 | `{1,2,3,5,6,7}` | `{4}` | −1 |
| `NLGF_P301S` | 3 | 6 | `{1,2,3,4,5,6}` | `{7}` | −1 |
| `NLGF_P301S` | 4 | 3 | `{2,3,4}` | `{1,5,6,7}` | −4 |
| **Total** |  |  |  |  | **−21** |

| Genotype | Observed / ideal AOIs |
|---|---:|
| `MAPTKI` | 20/28 |
| `P301S` | 20/28 |
| `NLGF_MAPTKI` | 28/28 |
| `NLGF_P301S` | 23/28 |
| **Total** | **91/112** |

Only `NLGF_MAPTKI` is fully sampled (28/28); `MAPTKI` `bio_rep = 1` is wholly
absent. Current loaders/descriptors exclude no AOIs: the gap predates this repo's live
GeoMx path. `storage/data/` contains only the processed `geomx.rds`, not the GeoMx
DCC/PKC/negative-probe/scan/mask inputs or missing-AOI records.

## snRNAseq — 16-unit structure, verified by contract

`R/core/io.R::load_snrnaseq()` loads the ~8.3G full object, subsets
`broad_annotations == "Microglia"` (~26k cells), and retains microglia + RNA counts.
Its hard contract requires:

- 4 genotypes × 4 batches = 16 populated (`genotype`, `batch`) units;
- exactly 16 `genotype_batch` IDs bijecting that fully crossed design; and
- no missing value in any design column.

Thus the report's 16/16 snRNAseq replicate units are verified-by-contract; a violated
cross or malformed identifier fails the build loudly.

## Known versus user-needed provenance

| Scope | Known in repo | Needed from user records |
|---|---|---|
| 24M bulk | 67-row manifest; first 16 rows = balanced Naoto 24M 2×2; exact raw labels, run order, canonical mapping, and fail-loud selection contract | Animal/aliquot IDs and lineage; 24M survival/attrition, litter, cage, sex, and batch records |
| Naoto 30M + Set6 | Exact unused labels, counts, row ranges, and acquisition-set membership | Why/disposition of unused arms; animal-level series confirmation; `MAPT-WT` disambiguation |
| Bulk assay method | TiO2 phospho acquisition/report identity documented in P7.2 | Spectronaut template and TiO2 enrichment/fractionation/acquisition method; see `.agent/p7_g2_bulk_assay_dossier.md` |
| GeoMx | Processed 91-AOI object; 112-cell ideal grid; exact 21-AOI deficit; retained design/QC metadata; no loader exclusion | DCC/PKC/negative-probe/scan/mask files; missing-AOI records; animal/aliquot linkage |
| snRNAseq | Microglia subset; fully crossed 4 × 4 design; 16 `genotype_batch` units verified by loader contract | Animal/aliquot linkage to GeoMx and bulk |
| Cross-assay | Gene-symbol harmonization only; no in-repo animal/aliquot crosswalk | Validated animal/aliquot crosswalk spanning snRNAseq ↔ GeoMx ↔ 24M bulk |

## Standing user-only record requests

1. Cross-assay animal/aliquot crosswalk spanning snRNAseq, GeoMx, and 24M bulk.
2. 24M animal survival/attrition, litter, cage, sex, and batch records.
3. GeoMx DCC/PKC/negative-probe/scan/mask inputs plus missing-AOI records.
4. Naoto 30M and Set6 arm disposition, animal-level series membership, and
   `MAPT-WT` disambiguation.
5. Spectronaut template plus TiO2 enrichment/fractionation/acquisition method; bulk-assay
   identity/request detail → `.agent/p7_g2_bulk_assay_dossier.md`.

## STANDING cross-modality paired-claim prohibition

Until a user-supplied cross-assay animal/aliquot crosswalk is validated:

- **PROHIBITED:** DIABLO, MOFA, mediation, within-animal correlation, or any other
  paired/animal-level cross-modality analysis or claim.
- **ALLOWED:** gene-symbol-level harmonization and unpaired modality-level
  meta-comparison, with each assay retaining its native replicate structure.

The current report uses only the allowed symbol-level comparison and modality-native
context; it is compliant with this gate.
