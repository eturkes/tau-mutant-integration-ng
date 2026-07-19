# P7.5 DAM-occupancy robustness results

## Decision

**FRAGILE.** The frozen primary E1 interaction remained positive in every variant, but
`primary_resolution=0.5` and `primary_resolution=0.6` crossed from the reference's supported
zero-null state to `fdr_zero > 0.05`. The pre-committed rule therefore rejects
**ROBUST-POSITIVE**. This is the result; positive direction in all variants and complete E2/E3
sign concordance do not override the frozen E1 decision rule.

## Frozen protocol

Executed `.agent/p7_dam_occupancy_prereg.md` frozen at `7755c9f`: reference plus 34
one-at-a-time variants; E1 beta-binomial probability-standardized, equal-present-batch DAM-fraction
interaction primary; E2 empirical-logit OLS plus batch-restricted Freedman-Lane permutation
(`seed=614`, `n_perm=9999`); E3 raw-proportion OLS; BH families over the five canonical contrasts;
margin `0.10`. Optional UCell re-scoring was not run because it is a declared optional secondary,
not part of the frozen 34-variant primary grid.

The work-unit cross-check was reconciled to the committed preregistration before execution. The
committed rule covers **all variants**, so any E1 `estimator_failed` would block
ROBUST-POSITIVE; sign tipping is `estimate <= 0`, including exact zero; margin reporting includes
both `fdr_minimum <= 0.05` and `abs(estimate) <= 0.10`. No estimator failed here, so the first
clarification does not alter this observed verdict.

## Reference anchor

| Check | Result |
|---|---:|
| Frozen substrate E1 interaction | `0.17411412574101759` (rounds to `0.1741141`) |
| Substrate vs frozen interaction-vector max abs diff | `0` |
| Cached membership path vs substrate max abs diff | `0` |
| Cached membership aggregation vs `microglia_state_substrate` | exact / `identical=TRUE` |
| Default re-annotation cell order | `identical=TRUE` |
| Default re-annotation labels | `identical=TRUE` |
| Re-annotation membership path vs substrate max abs diff | `0` |
| Grid | `35` variants = reference + `34` non-reference |
| Estimator attempts | `105` = `35 x 3` |

Reference E1 reproduced estimate `0.1741141`, 95% CI `[0.0954125, 0.2528158]`,
`fdr_zero=1.812987e-05`, `fdr_minimum=0.08116644`; E2 coefficient `0.9340592`, OLS 95% CI
`[0.1658742, 1.7022443]`, FDR `0.02806732`, permutation p `0.021`; E3 estimate
`0.1733368`, FDR `0.01854462`.

## Full variant x estimator table

`Status` order = E1/E2/E3. FDR0.10 = E1 minimum-effect family at margin `0.10`.

| Variant | Setting | Cells/units | E1 estimate [95% CI] | E1 FDR0 | E1 FDR0.10 | E2 coef [95% CI] | E2 FDR / perm p | E3 estimate [95% CI] | E3 FDR | Status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `reference` | reference labeling | 23160/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `resolution_0.2` | primary_resolution=0.2 | 23179/16 | 0.1802348 [0.0980381, 0.2624316] | 2.158e-05 | 0.089197 | 1.0538581 [0.2065449, 1.9011713] | 0.025322 / 0.0203 | 0.1775084 [0.0513144, 0.3037024] | 0.013936 | ok/ok/ok |
| `resolution_0.3` | primary_resolution=0.3 | 23168/16 | 0.1924237 [0.1055799, 0.2792674] | 1.759e-05 | 0.047633 | 1.0996647 [0.2040319, 1.9952975] | 0.026861 / 0.0208 | 0.1906944 [0.0511997, 0.3301890] | 0.016098 | ok/ok/ok |
| `resolution_0.5` | primary_resolution=0.5 | 23154/16 | 0.1090994 [-0.0017143, 0.2199130] | 0.089418 | 1.000000 | 0.8730018 [-0.3624002, 2.1084038] | 0.190235 / 0.1567 | 0.1154714 [-0.0558866, 0.2868295] | 0.269585 | ok/ok/ok |
| `resolution_0.6` | primary_resolution=0.6 | 23149/16 | 0.1213854 [0.0013640, 0.2414067] | 0.079088 | 1.000000 | 0.9240910 [-0.3816853, 2.2298673] | 0.186728 / 0.1614 | 0.1236647 [-0.0575092, 0.3048387] | 0.261605 | ok/ok/ok |
| `resolution_0.8` | primary_resolution=0.8 | 23084/16 | 0.1442624 [0.0172806, 0.2712443] | 0.043280 | 0.824145 | 1.1045101 [-0.2471159, 2.4561361] | 0.121966 / 0.1038 | 0.1546116 [-0.0354645, 0.3446876] | 0.164831 | ok/ok/ok |
| `id_floor_0.10` | id_floor=0.10 | 23160/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `id_floor_0.20` | id_floor=0.20 | 22798/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `mglike_floor_0.20` | mglike_floor=0.20 | 23915/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `mglike_floor_0.40` | mglike_floor=0.40 | 22798/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `no_prune` | no pruning | 26104/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `tol_0.05` | tol=0.05 | 23160/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `tol_0.15` | tol=0.15 | 23160/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `amb_floor_0.05` | amb_floor=0.05 | 23160/16 | 0.1777516 [0.0354968, 0.3200064] | 0.018378 | 0.473429 | 1.1782667 [-0.1585763, 2.5151097] | 0.096649 / 0.0795 | 0.2005870 [-0.0117124, 0.4128863] | 0.096545 | ok/ok/ok |
| `amb_floor_0.15` | amb_floor=0.15 | 23160/16 | 0.1741141 [0.0954125, 0.2528158] | 1.813e-05 | 0.081166 | 0.9340592 [0.1658742, 1.7022443] | 0.028067 / 0.0210 | 0.1733368 [0.0428491, 0.3038244] | 0.018545 | ok/ok/ok |
| `LOU_MAPTKI_batch01` | drop MAPTKI_batch01 | 22073/15 | 0.1634760 [0.0843092, 0.2426428] | 6.478e-05 | 0.145083 | 0.8417789 [-0.0007786, 1.6843364] | 0.062708 / 0.0487 | 0.1593056 [0.0150550, 0.3035561] | 0.042940 | ok/ok/ok |
| `LOU_MAPTKI_batch02` | drop MAPTKI_batch02 | 22276/15 | 0.1792319 [0.0966409, 0.2618229] | 2.633e-05 | 0.075093 | 1.0045297 [0.1480851, 1.8609744] | 0.033596 / 0.0206 | 0.1777503 [0.0294689, 0.3260317] | 0.030640 | ok/ok/ok |
| `LOU_MAPTKI_batch03` | drop MAPTKI_batch03 | 22079/15 | 0.1749919 [0.0910071, 0.2589766] | 5.538e-05 | 0.100127 | 0.9057259 [0.0332746, 1.7781771] | 0.054483 / 0.0452 | 0.1725753 [0.0238708, 0.3212797] | 0.035115 | ok/ok/ok |
| `LOU_MAPTKI_batch04` | drop MAPTKI_batch04 | 22127/15 | 0.1790988 [0.0946442, 0.2635533] | 4.041e-05 | 0.083008 | 0.9842024 [0.1182958, 1.8501091] | 0.038251 / 0.0397 | 0.1837159 [0.0374257, 0.3300060] | 0.025019 | ok/ok/ok |
| `LOU_NLGF_MAPTKI_batch01` | drop NLGF_MAPTKI_batch01 | 21747/15 | 0.1463134 [0.0660812, 0.2265455] | 0.000439 | 0.348509 | 0.7926626 [-0.0033352, 1.5886604] | 0.063450 / 0.0476 | 0.1467719 [0.0147640, 0.2787798] | 0.041804 | ok/ok/ok |
| `LOU_NLGF_MAPTKI_batch02` | drop NLGF_MAPTKI_batch02 | 21757/15 | 0.1395313 [0.0599635, 0.2190992] | 0.000735 | 0.461020 | 0.8271907 [-0.0038279, 1.6582094] | 0.063545 / 0.0530 | 0.1456647 [0.0151775, 0.2761519] | 0.041138 | ok/ok/ok |
| `LOU_NLGF_MAPTKI_batch03` | drop NLGF_MAPTKI_batch03 | 21195/15 | 0.2306348 [0.1668950, 0.2943745] | 1.654e-12 | 7.370e-05 | 1.1552160 [0.4908387, 1.8195932] | 0.004871 / 0.0018 | 0.2121185 [0.1020299, 0.3222072] | 0.002698 | ok/ok/ok |
| `LOU_NLGF_MAPTKI_batch04` | drop NLGF_MAPTKI_batch04 | 20602/15 | 0.1829263 [0.0932854, 0.2725671] | 7.931e-05 | 0.087260 | 0.9611675 [0.0884576, 1.8338775] | 0.043406 / 0.0413 | 0.1887919 [0.0455120, 0.3320718] | 0.020125 | ok/ok/ok |
| `LOU_NLGF_P301S_batch01` | drop NLGF_P301S_batch01 | 22927/15 | 0.1735197 [0.0881930, 0.2588464] | 8.408e-05 | 0.114083 | 0.9754972 [0.1065359, 1.8444584] | 0.040221 / 0.0271 | 0.1722522 [0.0235611, 0.3209433] | 0.035376 | ok/ok/ok |
| `LOU_NLGF_P301S_batch02` | drop NLGF_P301S_batch02 | 20521/15 | 0.1245674 [0.0594541, 0.1896807] | 0.000221 | 0.688067 | 0.7125811 [0.0489159, 1.3762463] | 0.047940 / 0.0440 | 0.1381491 [0.0203075, 0.2559907] | 0.033665 | ok/ok/ok |
| `LOU_NLGF_P301S_batch03` | drop NLGF_P301S_batch03 | 20936/15 | 0.2055501 [0.1326900, 0.2784102] | 4.017e-08 | 0.005651 | 1.0412773 [0.2105580, 1.8719966] | 0.025230 / 0.0168 | 0.2013659 [0.0713862, 0.3313456] | 0.009084 | ok/ok/ok |
| `LOU_NLGF_P301S_batch04` | drop NLGF_P301S_batch04 | 20458/15 | 0.1850716 [0.1018565, 0.2682867] | 1.633e-05 | 0.056379 | 1.0068814 [0.1517454, 1.8620173] | 0.033055 / 0.0315 | 0.1815798 [0.0343888, 0.3287709] | 0.027066 | ok/ok/ok |
| `LOU_P301S_batch01` | drop P301S_batch01 | 22208/15 | 0.1819643 [0.1019413, 0.2619874] | 1.040e-05 | 0.063255 | 1.0246134 [0.1808119, 1.8684149] | 0.028991 / 0.0188 | 0.1847859 [0.0390273, 0.3305444] | 0.023987 | ok/ok/ok |
| `LOU_P301S_batch02` | drop P301S_batch02 | 22222/15 | 0.1691552 [0.0853256, 0.2529849] | 9.571e-05 | 0.132383 | 0.8899201 [0.0218434, 1.7579967] | 0.057088 / 0.0442 | 0.1702347 [0.0217325, 0.3187369] | 0.036942 | ok/ok/ok |
| `LOU_P301S_batch03` | drop P301S_batch03 | 22598/15 | 0.1620481 [0.0827490, 0.2413471] | 7.746e-05 | 0.156414 | 0.7917872 [-0.0031731, 1.5867474] | 0.063404 / 0.0462 | 0.1618226 [0.0160981, 0.3075471] | 0.042010 | ok/ok/ok |
| `LOU_P301S_batch04` | drop P301S_batch04 | 21674/15 | 0.1825883 [0.1007731, 0.2644035] | 1.524e-05 | 0.064686 | 1.0299163 [0.1900173, 1.8698152] | 0.027786 / 0.0288 | 0.1765038 [0.0280107, 0.3249969] | 0.031763 | ok/ok/ok |
| `LOBO_batch01` | drop batch01 | 19475/12 | 0.1536186 [0.0633390, 0.2438981] | 0.001066 | 0.305502 | 0.8577956 [-0.1599299, 1.8755210] | 0.105969 / 0.0838 | 0.1506631 [-0.0241222, 0.3254485] | 0.099324 | ok/ok/ok |
| `LOBO_batch02` | drop batch02 | 17296/12 | 0.1203776 [0.0440833, 0.1966719] | 0.002481 | 0.899316 | 0.7075478 [-0.1731654, 1.5882609] | 0.121141 / 0.1231 | 0.1271756 [-0.0197430, 0.2740941] | 0.098112 | ok/ok/ok |
| `LOBO_batch03` | drop batch03 | 17328/12 | 0.2177949 [0.1488777, 0.2867122] | 7.334e-10 | 0.001010 | 1.0523863 [0.2465570, 1.8582156] | 0.023380 / 0.0175 | 0.2142382 [0.0797603, 0.3487161] | 0.010001 | ok/ok/ok |
| `LOBO_batch04` | drop batch04 | 15381/12 | 0.2028403 [0.1024055, 0.3032750] | 9.432e-05 | 0.063873 | 1.1185073 [0.0132750, 2.2237395] | 0.060057 / 0.0600 | 0.2012701 [0.0114526, 0.3910875] | 0.051201 | ok/ok/ok |

## E1 range and tipping set

- Empirical E1 interaction range: **[0.1090993551, 0.2306347619]**.
- Minimum: `resolution_0.5`; maximum: `LOU_NLGF_MAPTKI_batch03`.
- Direction tipping (`estimate <= 0`): **none**.
- Zero-null tipping (`fdr_zero > 0.05`):
  - `resolution_0.5`: estimate `0.1090993551`, `fdr_zero=0.0894183566`;
  - `resolution_0.6`: estimate `0.1213853525`, `fdr_zero=0.0790876621`.

## 0.10-margin behavior

E1 `fdr_minimum` range = **[7.369534e-05, 1.000000]**. The reference remains unresolved
(`0.08116644`). Four variants resolve `fdr_minimum <= 0.05`:

| Variant | E1 estimate | `fdr_minimum` |
|---|---:|---:|
| `resolution_0.3` | 0.1924237 | 0.04763314 |
| `LOU_NLGF_MAPTKI_batch03` | 0.2306348 | 7.369534e-05 |
| `LOU_NLGF_P301S_batch03` | 0.2055501 | 0.005650892 |
| `LOBO_batch03` | 0.2177949 | 0.001010092 |

Variants with `abs(E1 estimate) <= 0.10`: **none**; the empirical minimum is `0.1090993551`.
These descriptive margin results do not change the zero-null robustness verdict.

## Estimator failures

**0 / 105 estimator attempts failed.** Every variant was membership-estimable; all E1, E2, and
E3 fits returned `ok`. Failure inventory is empty.

## Cross-estimator concordance

- E2 sign versus E1: `35/35` concordant, `0` discordant.
- E3 sign versus E1: `35/35` concordant, `0` discordant.

This is sign concordance only. E2/E3 inferential support weakens at the same resolution variants,
but neither estimator defines the frozen verdict.

## Interpretation

The DAM-fraction interaction's estimated direction is stable across this frozen perturbation set,
and pruning/most annotation thresholds leave the reference estimate unchanged. However, the
report's central positive zero-null claim is **not robust under the pre-committed rule** because two
plausible clustering resolutions lose E1 BH-FDR support. The scientifically honest conclusion is
**FRAGILE**, specifically resolution-sensitive, with no estimator failure and no sign reversal.
