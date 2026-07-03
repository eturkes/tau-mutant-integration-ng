# Cross-Modality Narrative Figure Manifest

Scope: S1 audit for active plan "Cross-modality narrative figure pass".
Inputs read: `_crossmodality.qmd`, `_story.qmd`, `R/figures.R`, `_targets.R`,
`crossmodality_figures`, `story_figures`, rendered `_report/analysis-review.pdf`
pages 18-21 plus front story page. `tar_meta(crossmodality_figures, story_figures,
report)` = no target error/warning. Rscript emitted host `/etc/localtime` advisory
outside target metadata; not a figure defect.

Intent: preserve earned claims, replace dashboard-first integration figures with
named biological evidence plates. Reuse compact target data; add one axis-effect
spine in S2 before redrawing in S3.

## Current Figure Disposition

| fig_no | figure_id | current_target_slot | current_role | disposition | replacement_grammar | claim_boundary |
|---:|---|---|---|---|---|---|
| 5 | `fig-story-mechanism-crossmodality` | `story_figures$mechanism_crossmodality` (`pathway_axes`, `mechanism`, `clearance`) | Front synthesis: pathway-axis support + Myc/NF-kB/Gsk3b triage + clearance pairs. | REPLACE | 3-panel story mechanism plate: Myc TF effect forest; NF-kB/Gsk3b boundary mini-forest; focused Apoe-Trem2 / clearance effect strip from measured modalities. Avoid broad status matrices. | RNA supports Myc interaction; NF-kB attenuation discordant / unsupported; bulk phospho does not recover Gsk3b; Apoe-Trem2 support is focused, not full CCC. |
| 31 | `fig-crossmodality-four-modality-counts` | `crossmodality_figures$four_modality_counts` | Global assay-family FDR count heatmap by contrast. | DEMOTE | Supplemental/audit count strip only if retained; remove from main biological flow. | Counts show assay burden and amyloid-dominant scale; they do not establish axis biology or broad "all modalities support" claims. |
| 32 | `fig-crossmodality-four-modality-pathways` | `crossmodality_figures$four_modality_pathways` | Pathway-axis x contrast support matrix. | REPLACE | Split into named-axis evidence subpanels: amyloid->DAM/AP, synaptic suppression, clearance context, NF-kB boundary. Use effect/direction dot or lollipop rows by modality, not a single all-axis heatmap. | Direction can be mixed; NF-kB axis support is context, not a supported attenuation mechanism. |
| 33 | `fig-crossmodality-four-modality-symbols` | `crossmodality_figures$four_modality_symbols` | Selected symbol x modality x contrast status heatmaps. | REPLACE | Axis-symbol effect plate: selected genes/pathways as rows; modality columns carry signed effect, FDR/support glyph, and measured/unmeasured state. Facet by biological axis, not contrast dashboard. | Selected symbols illustrate named axes only; unmeasured/absent stays explicit; bulk rows are hippocampal, not microglia-sorted. |
| 34 | `fig-crossmodality-geomx-counts` | `crossmodality_figures$geomx_counts` | GeoMx up/down FDR count stems by focal contrast. | DEMOTE | Convert to a small boundary/count inset or fold into GeoMx volcano captions; no standalone main figure. | Spatial DE burden is largest for amyloid, but count magnitude alone is not a cell-state or pathway claim. |
| 35 | `fig-crossmodality-geomx-volcano` | `crossmodality_figures$geomx_volcano` (`points`, `labels`, `counts`) | Conventional GeoMx DE volcanoes. | KEEP-RESHAPE | Use as GeoMx row in amyloid-response concordance plate; label DAM/AP/synaptic/clearance symbols selected by S2 rules; optionally limit facets to amyloid-on-MAPTKI/P301S plus interaction boundary. | GeoMx supports dense amyloid spatial DE; SpatialDecon abundance remains blocked; interaction support remains bounded. |
| 36 | `fig-crossmodality-geomx-sensitivity` | `crossmodality_figures$geomx_sensitivity` | GeoMx replicate-treatment sensitivity audit. | DEMOTE | Compact robustness/boundary inset after primary spatial evidence. | Primary GeoMx model is load-bearing; unblocked and bio-unit-collapsed fits expose sensitivity, not new biology. |
| 37 | `fig-crossmodality-bulk-counts` | `crossmodality_figures$bulk_counts` | Bulk proteome/phospho FDR count stems. | DEMOTE | Replace standalone counts with selected effect rows inside amyloid-response and synaptic-clearance plates. | Bulk corroborates context but is not microglia-sorted; count totals cannot localize effects to microglia. |
| 38 | `fig-crossmodality-phospho-correction` | `crossmodality_figures$phospho_raw_corrected` (`points`, `labels`, `counts`) | Raw-vs-parent-corrected phosphosite effect geometry. | KEEP-RESHAPE | Keep as a boundary/mechanism panel; label selected phosphosites/kinase-relevant rows and tie to interaction/Gsk3b non-recovery. | Parent correction changes geometry; 24M hippocampal phospho covers the layer but does not recover Gsk3b interaction or tau-in-NLGF support. |

## S2/S3 Target-Slot Handoff

Add one compact spine first, then draw plates from it.

| planned_slot | producer_target | contents | needed_by |
|---|---|---|---|
| `axis_effect_spine` | `crossmodality_figures` | Long table: `axis`, `feature_id`, `feature_label`, `feature_type` (`symbol`/`pathway`/`pair`/`boundary`), `modality_class`, `contrast`, `effect`, `effect_scale`, `fdr`, `support_status`, `direction`, `measured_state`, `source_slot`, `selection_rank`. Complete keys for selected axes x modalities x focal contrasts. | All replacement plates; S2 tests complete keys, finite measured effects, deterministic order, explicit blocked/unmeasured states. |
| `amyloid_response_plate` | `crossmodality_figures` | Subset of `axis_effect_spine` for DAM/AP amyloid response plus `geomx_volcano` labels/counts. | Replace figures 31-35 main integration opening. |
| `synaptic_clearance_plate` | `crossmodality_figures` | Synaptic-loss and clearance/DAM features, including Apoe-Trem2 measured-pair support and unearned pairs. | Replace figure 33 status matrix and part of story Figure 5. |
| `interaction_boundary_plate` | `crossmodality_figures` | Interaction rows across microglia composition, GeoMx, proteome, raw/corrected phospho, Myc/NF-kB/Gsk3b boundary fields. | Replace count/status treatment of interaction; carry negative evidence as effects/CIs/status, not boxes alone. |
| `mechanism_crossmodality` | `story_figures` | Redrawn front story plate sourced from `axis_effect_spine`, existing `mechanism`, and `clearance` rows. | Replace Figure 5 while preserving target slot name. |

Selection rules for S2:
- Axes fixed: `DAM`, `antigen_presentation`, `synaptic`, `clearance`, `interaction_boundary`, `mechanism_boundary`.
- Focal contrasts fixed: `nlgf_in_maptki`, `nlgf_in_p301s`, `tau_in_nlgf`, `interaction`.
- Feature priority: predeclared project anchors first (`Apoe`, `Trem2`, `Cst7`, `Spp1`, `Cd74`, `Syn1`, `Syp`, `Snap25`, `Gsk3b`, `Myc`, NF-kB family rows), then top deterministic ranked rows from existing compact pathway/symbol tables.
- Encoding contract: `measured_state` in `measured`, `not_observed`, `blocked`, `not_applicable`; unmeasured/blocked never coerced to non-significant.

## Acceptance Check

S1 acceptance satisfied when this file exists and roadmap marks S1 done:
- Exact figure IDs: Figure 5 plus Figures 31-38 listed above.
- Exact current target slots: listed in `current_target_slot`.
- Biological role, disposition, replacement grammar, claim/boundary: listed per row.
- S2 target handoff: exact planned slots and selection/encoding contract listed.
