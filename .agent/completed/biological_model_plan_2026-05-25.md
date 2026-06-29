# Active plan: integrated biological model synthesis (Phase H)

- **Created:** 2026-05-24 (immediate successor to
  `mechanism_layer_plan_2026-05-24.md`).
- **Last updated:** 2026-05-25 (H5 DONE; upfront-mirror snippet
  inserted as section 1.1 of analysis.Rmd; plan moves to
  `storage/notes/completed/biological_model_plan_2026-05-25.md` at
  this commit; Phase H COMPLETE).
- **State:** Phase H COMPLETE. H1 DONE 2026-05-25. H2 DONE 2026-05-25
  (75 rows, Strong=24 / Moderate=39 / Suggestive=12; net support
  T-Synergy 28 > Hyp-3B 27 > T-Inflammation 26 > Hyp-1B 24). H3 DONE
  2026-05-25 (per-axis contest verdicts Hyp-1B / Hyp-2B / Hyp-3B win their
  contests by net-support margins 6 / 12 / 53 respectively, no
  tie-break rule invoked; cross-axis theme ranking by net_support
  is T-Synergy 28 > T-Inflammation 26 > T-Compartment-suppression 17
  > T-Tau-attenuates 5, with Hyp-0 stand-alone 7). H4 DONE 2026-05-25
  (chapter renders as section 17 with 5 subsections; "Session info"
  cleanly renumbers from section 17 to section 18 under the parent's
  number_sections: true; single-axis anchoring + Hyp-0 evenness named
  in 17.3.2; integrated reading hedged as "currently best-supported",
  not authoritative). H5 DONE 2026-05-25 (upfront-mirror inserted as
  section 1.1; ~205 words across 3 paragraphs; zero errors / zero
  warnings in re-knitted analysis.html; see Outcome summary at the
  bottom of this file for the per-axis verdict + ledger / adjudication
  / chapter / mirror manifest + next-plan pointer).
- **Predecessor:** `storage/notes/completed/mechanism_layer_plan_2026-05-24.md`
  (Phase G COMPLETE; the mechanism arc closed with per-axis TF + kinase
  + LR verdicts at sections 14.3 / 15.3 / 16.3 of the knitted
  `analysis.html`).
- **Goal:** synthesise the project's three-axis verdict (D2 cross-
  modality leader board) and three mechanism-layer verdicts (E5 TF /
  F5 kinase / G5 LR) into an explicit integrated biological model
  written into the R Markdown. The synthesis takes a **structured
  hypothesis-comparison** form: a confidence-graded claims ledger
  adjudicates a set of competing biological models on the project's
  own evidence; the late chapter renders the ledger + the adjudication
  + the working integrated model; a short upfront mirror in
  `analysis.Rmd` echoes the headline conclusion at the top of the
  document. The synthesis adds no new computation — every claim is
  derived from existing caches / TSVs / verdicts / completed-plan
  outcome summaries.

## Why this plan exists

`analysis.html` now carries 16 sections of completed analysis plus
three per-axis mechanism-layer verdicts, but there is no integrated
reading anywhere in the document. The closest substitutes are
`storage/notes/completed/mechanism_layer_plan_2026-05-24.md` §
"Final cross-axis verdict" (an internal-only synthesis at line 3421)
and the individual verdict TSVs (axis-by-axis prose paragraphs in
`evidence_summary`, not cross-axis). A reader of `analysis.html`
itself sees the per-section findings but never reads "what is the
biology", "which competing hypotheses survive these data", or "what
are the testable predictions of the favoured model".

The synthesis chapter and upfront mirror close that gap. The
hypothesis-comparison stance is deliberate: a strictly narrative
synthesis would be vulnerable to LLM-style training-data bias toward
canonical AD-microglia stories (classical complement-mediated
synaptic pruning; uniform amyloid-tau amplification); the structured
claims-ledger + per-model-adjudication form keeps the synthesis
data-grounded and the competing models explicit. The confidence-
grading column makes the synthesis honestly stratify weak from strong
claims.

## Locked decisions (carried forward from `mechanism_layer_plan_2026-05-24` + H0 user direction + H1)

| Decision | Choice |
|---|---|
| Synthesis venue | **New `rmd/16_biological_model.Rmd` rendering as section 17** + **short upfront "Headline integrated model" snippet** inserted near the top of `analysis.Rmd` after the scientific-question subsection and before the data section (locked at H0 user direction 2026-05-24; selected from the four-option AskUserQuestion: late chapter + upfront mirror). Late chapter is the primary venue (the structured ledger + model adjudication live there); upfront snippet is 1-2 paragraphs naming the contests and pointing forward to section 17 for the ledger. |
| Synthesis form | **Confidence-graded claims ledger only** (locked at H0 user direction 2026-05-24; selected from the four-option AskUserQuestion: "Confidence-graded claims ledger only" over "Standard / Lean / Deep with schematic figures"). The synthesis chapter is structured around a per-claim ledger; per-axis narrative prose is intentionally avoided in favour of the structured form. Schematic figures are out of scope for this plan. |
| Synthesis stance | **Hypothesis comparison** (locked at H0 user direction 2026-05-24; selected from the three-option AskUserQuestion: "Hypothesis comparison" over "Working model with confidence levels / Authoritative model"). The synthesis frames the integrated reading as a set of competing biological models adjudicated on the project's evidence rather than as a single working model. The favoured model emerges per axis from the adjudication, never asserted up-front. |
| New computation | **None** (locked at H0 by construction). Every claim in the ledger is derived from existing verdict TSVs, the pathway leader board, integration tables, completed-plan outcome summaries, and the prose narratives in the existing Rmd sections. No new caches, no new DE runs, no new tool calls. The synthesis is "structure what is already there", not "compute more". |
| Competing-model set | **Hybrid: per-axis contests + cross-axis themes overlay (H1 decision gate locked 2026-05-25 via user steer "Both of these sound interesting", refined into a hybrid via follow-up AskUserQuestion).** Eleven entities total: (a) six per-axis binary-contest models = Hyp-1A/Hyp-1B at axis 1 amyloid_activation; Hyp-2A/Hyp-2B at axis 2 synaptic_suppression; Hyp-3A/Hyp-3B at axis 3 interaction_metabolic; (b) Hyp-0 Cdk5 cross-axis integrator (stand-alone, no contest); (c) four cross-axis themes T-Inflammation, T-Compartment-suppression, T-Tau-attenuates, T-Synergy (stand-alone, no head-to-head contest with each other). Per-axis contests adjudicated binary (winner per axis); Hyp-0 + themes adjudicated as independent stand-alone evidence counts. Full model definitions and theme statements in the H1 step body (DONE section below). Claims can support / contradict any of the 11 entities in their ledger row; a single claim commonly supports both a per-axis model (e.g. Hyp-1B) and the theme that subsumes it (e.g. T-Tau-attenuates). |
| Claims-ledger schema | **14 columns with full audit trail (H1 decision gate locked 2026-05-25; user picked the recommended default over the lean 8-column alternative).** Columns: `claim_id`, `axis`, `layer`, `claim` (one-sentence text), `direction` (+ / − / mixed), `effect_size` (concrete number string), `primary_evidence_source` (TSV path or section ref), `corroborating_evidence` (semicolon-joined references), `n_replicates_in_modalities`, `n_replicates_in_layers`, `confidence_grade` (Strong / Moderate / Suggestive), `supports_models` (semicolon-joined model IDs from the 11-entity set), `contradicts_models` (semicolon-joined model IDs from the 11-entity set), `notes` (free-text caveats including Interpretation A entanglement, sign-reversal, n_cells_used = 1, etc.). Saved as `storage/results/biological_model_claims_ledger.tsv`. The audit trail is what lets H4 defend a Suggestive grade when a reader clicks the row to ask "why suggestive?" — the schema is designed so adjudication in H3 is a structured groupby on `supports_models` / `contradicts_models` rather than free-text parsing. |
| Confidence-grade definitions | **Three-tier Strong / Moderate / Suggestive (H1 decision gate locked 2026-05-25; user picked the recommended default over the binary Strong / Weak alternative).** **Strong** = (a) effect reproduces across ≥ 2 modalities with consistent sign AND (b) reaches FDR < 0.10 in primary inference tool AND (c) has cross-axis or cross-layer corroboration named explicitly in the corroborating-evidence column. Example: Gsk3b at interaction (padj 0.00152, kinase layer, corroborates Myc TF suppression and axon-guidance LR rewiring at the same axis). **Moderate** = effect is significant in primary inference tool (FDR < 0.10) but lacks ≥ 2-modality replication OR lacks cross-layer corroboration. Example: Apoe_Trem2 at axis-2 rank 9 in LIANA+ without cross-tool ranking at the same axis-relevant contrasts (the G3-vs-G4 complementarity property). **Suggestive** = effect is in top-N of primary inference tool but does not reach significance OR has known structural caveats (e.g. Interpretation A entanglement at axes 1 + 2 for TF and kinase layers; single-cell n_cells_used = 1 LR cells). Example: most TF / kinase top-5 at axis 2 inherit from axis-1 ranking under Interpretation A entanglement. |

## Open questions (defer to the session that needs them)

| Question | Default proposal | Decided in |
|---|---|---|
| Ledger source-priority order | Verdict TSVs first (E5/F5/G5); then pathway leader board (D2); then cross-modality integration (D1); then per-state NEBULA; then completed-plan outcome summaries. Source most-distilled to least. | H2 (executes the priority). |
| How to handle Interpretation A entanglement at axes 1 + 2 | Note explicitly in the ledger `notes` column; downgrade affected claims to Suggestive when the entanglement is the load-bearing replication argument; treat as Moderate when in-universe-target-count differential adds independent evidence. | H2 (per-claim judgment). |
| How to handle the cross-tool vs axis-restricted Hyp-2A nuance | The `Pros1_Mertk` LR is absent from the axis-2 restricted leader board (Hyp-2A appears contradicted) but appears as a strong cross-tool LR (leader_score 15.4, sign +, several microglia↔microglia combinations) in `ccc_lr_cross_tool_leaderboard.tsv` (Hyp-2A partially supported when scope widens). Emit two ledger rows: one row for the axis-restricted view (contradicts Hyp-2A, supports Hyp-2B, supports T-Compartment-suppression-as-non-complement-route); one row for the cross-tool view (supports Hyp-2A's classical-complement axis specifically as Pros1-Mertk efferocytosis; notes column flags scope difference). The H4 17.3 adjudication paragraph names the scope ambiguity honestly rather than collapsing it. | H2 (per-claim judgment, with verification subagent finding 2026-05-25). |
| Adjudication arithmetic | Per-model: count rows in the ledger where `supports_models` includes the model; count rows where `contradicts_models` includes it; emit `support_count`, `contradict_count`, `net_support` ( = support − contradict), and `support_to_contradict_ratio`. Per-(model, axis): same arithmetic restricted to the axis. The favoured model per contest is the higher-ranked by net_support; ties broken by `n_strong_supports`. Hybrid scheme: per-axis contests get a binary favoured-model verdict; Hyp-0 + the four cross-axis themes get stand-alone support-count summaries with no head-to-head winner. Theme adjudication uses the same arithmetic as per-axis but reports `support_count`, `contradict_count`, `net_support`, and a confidence-distribution histogram (n Strong / Moderate / Suggestive supports). | H3. |
| Upfront mirror length | 2-3 short paragraphs: (1) name the three per-axis contests; (2) state the favoured-model verdict per axis in one clause each + name the cross-axis-theme leaders; (3) point to section 17 for the ledger + adjudication + integrated model statement. ~150 words. | H5. |

## Execution model

- A fresh session must read `CLAUDE.md` first, then this file in full.
- A fresh session executes the NEXT step whose status is `TODO`. Steps
  are sized to fit comfortably in one session each.
- After completing a step, the session must:
  1. Update that step's status from `TODO` to `DONE <YYYY-MM-DD>` with a
     multi-paragraph completion note (what was built, key file paths,
     biology highlights, plan-spec deviations explicitly acknowledged).
  2. `chown rstudio:rstudio` any new files (agent runs as root).
  3. Re-knit `analysis.Rmd` and verify zero `class="error"` and zero
     `class="warning"` in the rendered HTML.
  4. Commit locally per CLAUDE.md (imperative subject under 70 chars;
     co-author trailer).
  5. Either continue to the next `TODO` step in the same session if the
     context budget permits, or end cleanly.
- Decision-gate step H1 requires user confirmation via
  `AskUserQuestion`. Wait for confirmation before proceeding to H2.

## Phase H: integrated biological model synthesis (sessions H1-H5)

### Session H1: model framework + ledger schema decision gate [DONE 2026-05-25]

**Completion note (2026-05-25).** H1 ran as a pure decision-gate
session per the plan spec ("no new R code, no new TSV, no Rmd
modification"). The session unfolded in three movements: (i) a
verification subagent (`general-purpose`) confirmed that all the
quantitative anchors of the default model proposals reproduce
exactly in the existing TSVs — Rela sign-reversal (+2.565, +1.809 at
NLGF arms, −1.356 at interaction); Spi1 and Nfkb1 axis-1 activations
(ranks 2 and 5); Apoe_Trem2 axis-2 leader_rank 9 + 10; App_Cd74
axis-2 rank 12 + 20; classical complement absent from axis-2
restricted leader board; 49/100 top-100 axis-2 microglia-involving;
Gsk3b padj = 0.00152 at interaction with per-contrast values
matching the proposal exactly (amyloid_alone −1.660, tau_alone
−1.469, interaction +4.323); Myc/Creb1/Tp53/Jun all suppressed at
interaction (ranks 2-5, signs −,−,−,−); Cdk5 at top-5 of all three
axes (ranks 2, 2, 4) with substrate footprints 81 / 263 / 153 in
`n_targets_in_axis_universe`. (ii) A material nuance surfaced during
verification: `Pros1_Mertk` is absent from the axis-2 restricted
leader board (supporting the Hyp-2A contradiction the proposal
predicted) but IS strong in `ccc_lr_cross_tool_leaderboard.tsv`
(leader_score 15.4, sign +) — so Hyp-2A's contradiction is precise only
when scoped to the axis-restricted view. This nuance was added to
the Open-questions table as an explicit H2 per-claim judgment with a
two-row resolution protocol (axis-restricted row vs cross-tool row).
(iii) An initial three-question AskUserQuestion presented the
default + alternative for each H1 lock; the user picked the
recommended defaults for the schema (14 columns) and grades (3-tier)
but answered "Both of these sound interesting" for the model-set
question. A follow-up single-question AskUserQuestion presented
three merge variants — hybrid (recommended), flat-themes-only, and
per-axis-only with prose rollup — and the user picked the hybrid.

The hybrid scheme is the structural payload of H1: the H2 ledger
will tag each claim with support / contradict membership in any of
**eleven** named entities rather than the six originally
proposed. Six are the per-axis binary contests (Hyp-1A/Hyp-1B at axis 1
amyloid_activation; Hyp-2A/Hyp-2B at axis 2 synaptic_suppression;
Hyp-3A/Hyp-3B at axis 3 interaction_metabolic). One is the cross-axis
integrator slot (Hyp-0 Cdk5). Four are cross-axis themes that span
multiple axes and are intentionally allowed to subsume per-axis
models (e.g. T-Tau-attenuates subsumes Hyp-1B; T-Synergy subsumes
Hyp-3B; T-Compartment-suppression subsumes Hyp-2A; T-Inflammation
subsumes both Hyp-1A and Hyp-2A jointly). Per-axis contests retain their
binary head-to-head adjudication; Hyp-0 and the four themes are
adjudicated as stand-alone evidence counts (no winner). H3 will emit
two adjudication outputs: a per-axis-contest verdict table and a
cross-axis-theme support tally. The expectation is that a single
claim row commonly populates both a per-axis model column AND the
theme column it inherits — e.g. the Rela sign-reversal claim
supports Hyp-1B and T-Tau-attenuates simultaneously.

The full 11-entity definitions are spelled out in the next two
subsections below; H2 should consult them before assigning
`supports_models` / `contradicts_models` cells. The H2 claim-
categories handoff sketch is the third subsection. All H1 deliverables
are plan-edit only; no R or Rmd files were touched, so the Phase G
analysis.html (Phase G COMPLETE per the predecessor plan's outcome
summary) remains the current rendered state.

#### H1 deliverable 1: full 11-entity model definitions (LOCKED 2026-05-25)

The H2 ledger's `supports_models` / `contradicts_models` columns
accept any subset of these eleven entity IDs as semicolon-joined
values. Definitions are grouped by structural type.

**Per-axis binary contests (6 models, 3 contests).** Per-axis
contests produce a binary favoured-model verdict in H3.

- **Contest 1 (axis: amyloid_activation)** — adjudicates how
  amyloid reshapes microglia and whether tau modifies the amyloid
  program.
  - **Hyp-1A — Tau-independent amyloid program.** Amyloid drives a
    canonical DAM-state inflammatory program; the program is
    qualitatively identical in MAPTKI and P301S backgrounds; tau
    contributes no inflammatory modulator. The interaction contrast
    at axis 1 is zero or near-zero. Predicts: matched-direction
    matched-magnitude TF / kinase / LR / pathway signals at
    `nlgf_in_maptki` and `nlgf_in_p301s` with no interaction-contrast
    signal.
  - **Hyp-1B — Tau attenuates amyloid-driven NF-κB.** Amyloid drives
    NF-κB-led inflammation in MAPTKI; tau on the amyloid background
    SUPPRESSES this program. Rela sign-reversal is the load-bearing
    signature: positive at `nlgf_in_maptki` / `nlgf_in_p301s`,
    negative at `interaction`. The interaction contrast at axis 1 is
    biologically informative. Predicts: Rela score sign-flip; Nfkb1
    score-attenuation at interaction; coordinated Spi1 / Nfkb1 / Rela
    asymmetric reading at the two NLGF arms.

- **Contest 2 (axis: synaptic_suppression)** — adjudicates whether
  NLGF-driven synaptic suppression is mediated by classical
  complement or by TREM2/APP/synaptic-adhesion routes.
  - **Hyp-2A — Classical complement-pruning hypothesis.** NLGF-driven
    synaptic suppression is mediated by canonical C1q / C3 /
    Mertk-Pros1 / Tyro3-Axl complement-tagging + receptor-mediated
    engulfment. Predicts these LR pairs at the top of axis-2 cell-
    cell communication rankings.
  - **Hyp-2B — TREM2 / APP-fragment-mediated clearance.** NLGF-driven
    synaptic suppression is mediated by TREM2-mediated clearance +
    APP-fragment uptake + synaptic-adhesion modulation (Apoe_Trem2,
    App_Cd74, Cd200_Cd200r1, Vcan_Cd44, Ncam1 ligands). Predicts
    these pairs at the top of axis-2 rankings instead of classical
    complement.

- **Contest 3 (axis: interaction_metabolic)** — adjudicates whether
  the tau × amyloid interaction produces a qualitatively new
  mechanism or is arithmetically the difference of single-insult
  arms.
  - **Hyp-3A — No interaction-specific mechanism.** The interaction
    contrast is the arithmetic difference of single-insult arms; no
    qualitatively new mechanism emerges at the tau × amyloid
    interface. Mechanism-layer signals at the interaction are
    smaller-than-amyloid-alone in magnitude or absent.
  - **Hyp-3B — Distinct synergy mechanism (Gsk3b / Myc / axon-guidance).**
    At the tau × amyloid interface a qualitatively new mechanism
    emerges: Gsk3b activation (the strongest single-kinase signal of
    the project, padj 0.00152, supported by Mapk14 / Cdk5 / Cdk1
    canonical tau-kinase activation), Myc / Creb1 / Tp53 / Jun
    coordinated biosynthetic suppression at the TF layer, and
    axon-guidance / synaptic-adhesion intercellular rewiring at the
    LR layer. The mechanism is NOT inferrable from either single-
    insult arm because all three layers carry interaction-contrast-
    specific sign or magnitude patterns.

**Cross-axis integrator slot (1 model, 0 contests).** Adjudicated
as a stand-alone evidence count.

- **Hyp-0 — Cdk5 cross-axis integrator.** Cdk5 is in the top 5 of the
  kinase layer at all three axes (ranks 2 / 2 / 4) and its KSN
  substrate footprint scales with the axis universe size (81 sites
  at amyloid, 263 at synaptic, 153 at interaction). Predicts:
  Cdk5-substrate enrichment in claims at every axis;
  CDK5R1/CDK5R2 expression continuity across genotypes; cross-
  axis TF or pathway concordance traceable to Cdk5-phosphorylation
  targets.

**Cross-axis themes (4 themes, 0 contests).** Adjudicated as
stand-alone evidence counts, in parallel; themes do not compete
head-to-head with each other. Themes intentionally overlap with
per-axis models — a claim row that supports a per-axis model
typically also supports the theme it subsumes, and the ledger
documents both.

- **T-Inflammation — Additive DAM-amplification across axes.** A
  unified inflammatory program drives both axes 1 and 2; amyloid
  and tau act additively on a canonical DAM inflammatory state;
  the inflammatory program is qualitatively the same regardless of
  tau background. Predicts: axis-1 TF / kinase activations replicate
  at axis-2 contrasts; no qualitative interaction-specific mechanism
  at axis 3. Theme is contradicted by Hyp-1B / Hyp-3B evidence (sign-
  reversal, novel synergy mechanism). Subsumes Hyp-1A and Hyp-2A jointly.
- **T-Compartment-suppression — Synaptic compartment is suppressed
  via classical complement.** NLGF suppresses the synaptic
  compartment via canonical complement-mediated pruning. Predicts:
  axis-2 evidence shows synaptic GO term suppression AND complement
  LR pairs at the top of axis-2 rankings. Theme subsumes Hyp-2A.
- **T-Tau-attenuates — Tau modifies amyloid rather than amplifying
  it.** Tau modifies the amyloid program rather than amplifying it;
  specifically, tau suppresses amyloid-driven NF-κB at the
  interaction (Rela sign-reversal). Predicts: positive Rela at NLGF
  arms, negative Rela at interaction; coordinated NF-κB-family
  asymmetric reading across contrasts. Theme subsumes Hyp-1B.
- **T-Synergy — Qualitatively new mechanism at the tau × amyloid
  interface.** A qualitatively new mechanism emerges at the
  interaction that is not present in either single-insult arm.
  Predicts: (a) Gsk3b activation specifically at interaction; (b)
  coordinated TF suppression (Myc/Creb1/Tp53/Jun) at interaction;
  (c) interaction-specific LR rewiring (axon-guidance / Gpc3_Unc5c /
  Sema6a_Plxna4 sign-flip). Theme subsumes Hyp-3B.

#### H1 deliverable 2: H2 claim-categories handoff sketch

H2 should target **40–80 ledger rows total**, distributed roughly
as follows. Each entry below names: claim category, primary
evidence source, expected confidence range, and the entities the
category most often supports (S) / contradicts (C). The sketch is
non-exhaustive; H2 finalises by reading the verdict TSVs row-by-row
and emitting one ledger row per atomic claim.

- **Axis 1 / TF layer** (~7–10 rows): Spi1 activation; Nfkb1
  activation; Sp3 activation; Rreb1 repression; A0A979HLR9
  activation; Rela activation at NLGF arms; Rela suppression at
  interaction (sign-reversal, separate row). Mostly Moderate;
  Rela sign-reversal claim is Strong (cross-axis corroboration
  via the kinase-layer interaction signature). S: Hyp-1A or Hyp-1B (per
  row), T-Inflammation, T-Tau-attenuates (for sign-reversal row).
  Sources: `tf_activity_axis_restricted.tsv`,
  `tf_activity_unified_leaderboard.tsv`.
- **Axis 1 / kinase layer** (~5–7 rows): Cdk2 activation; Cdk5
  activation (p301s-biased — also feeds Hyp-0); Mapk8 activation
  (p301s-biased); Csnk1e activation; Camk2g repression. Mostly
  Moderate; Cdk5 contributes to Hyp-0 Strong. Sources:
  `kinase_activity_axis_restricted.tsv`,
  `kinase_activity_per_contrast.tsv`.
- **Axis 1 / LR layer** (~8–12 rows): top-5 axon-guidance bias
  (Adam11_Itga4 / Gas6_Axl / L1cam_Egfr / Omg_Rtn4r / Rspo3_Lgr4);
  Cd200_Cd200r1; Apoe_Trem2 (rank 9); App_Cd74 (rank 11);
  Tgfb1_Sdc2 (rank 34); 55/100 microglia-involving claim;
  Pros1_Mertk cross-tool row (per the H2 nuance protocol). Mostly
  Moderate; Apoe_Trem2 at axis 1 contradicts Hyp-2A's predicted-
  location framing for axis 1 (Hyp-2A is an axis-2 model, so this is
  a context note). Sources: `ccc_lr_axis_restricted.tsv`,
  `ccc_lr_cross_tool_leaderboard.tsv`,
  `ccc_multinichenet_top100_*`.
- **Axis 1 / pathway layer** (~5–7 rows): MG-M2 module
  (cross-modality consistent-sign = 4 at nlgf_in_p301s);
  GOBP_ADAPTIVE_IMMUNE_RESPONSE; MHC II antigen-presentation
  family; custom DAM_up / WAM / HAM_disease / AD2. Mostly Strong
  (cross-modality replication built in). S: T-Inflammation, Hyp-1A.
  Sources: `pathway_survey_unified_leaderboard.tsv`,
  `hdwgcna_module_de.tsv`, `hdwgcna_module_enrichment.tsv`.
- **Axis 2 / TF layer** (~3–5 rows): the honest non-finding
  claim (no axis-2-restricted TF passes FDR with
  cross-modality replication); Interpretation A entanglement
  caveat; in-universe target-count differential. Mostly
  Suggestive; the non-finding is itself a claim that supports
  M-shared (Hyp-1A/Hyp-2A) by absence. Sources: verdict TSV
  `tf_activity_verdict.tsv`.
- **Axis 2 / kinase layer** (~3–5 rows): honest non-finding;
  Interpretation A entanglement; Cdk5 KSN footprint expansion
  (81 → 263 sites — the Hyp-0-supporting row). Mostly Suggestive;
  Cdk5 KSN-footprint row is Moderate-to-Strong (cross-axis,
  per Hyp-0). Sources: `kinase_activity_verdict.tsv`,
  `kinase_activity_axis_restricted.tsv`.
- **Axis 2 / LR layer** (~10–14 rows): Apoe_Trem2 (MG_DAM →
  MG_DAM rank 9); Apoe_Trem2 (Astrocyte → MG_DAM rank 10);
  App_Cd74 (Oligodendrocyte → Microglia_proliferative rank 12);
  App_Cd74 (Oligodendrocyte → MG_DAM rank 20); 49/100 top-100
  microglia-involving; 6,459 axis-2 cells (largest axis);
  C1qb_Lrp1 (rank 775, the first complement pair); absence of
  C3 / Mertk-Pros1 / Tyro3-Axl in top 800 (axis-restricted);
  Pros1_Mertk cross-tool row (the H2 nuance protocol). Multiple
  Strong rows here. S: Hyp-2B, T-Compartment-suppression;
  contradicts Hyp-2A on axis-restricted scope. Sources:
  `ccc_lr_axis_restricted.tsv`, `ccc_lr_cross_tool_leaderboard.tsv`.
- **Axis 2 / pathway layer** (~3–5 rows): six synaptic GO CC
  terms at leader_score 11.6–15.4 (consistent-sign negative at
  NLGF contrasts); the synaptic-suppression "what" claim that
  every contest interprets differently. Strong. S:
  T-Compartment-suppression, Hyp-2A AND Hyp-2B (the synaptic
  suppression itself is shared; only the mechanism differs).
- **Axis 3 / TF layer** (~5–7 rows): Myc suppression; Creb1
  suppression; Tp53 suppression; Jun suppression; A0A979HLR9
  interaction-context negative (sign-reversal vs axis 1).
  Mostly Strong (cross-layer corroboration via Gsk3b kinase
  signature). S: Hyp-3B, T-Synergy.
- **Axis 3 / kinase layer** (~5–7 rows): Gsk3b activation
  (THE strongest single-kinase result; padj 0.00152; per-
  contrast amyloid_alone −1.660, tau_alone −1.469, interaction
  +4.323; SYNERGY IS the activation); Mapk14 activation; Cdk5
  activation (cross-axis, feeds Hyp-0); Cdk1 activation; Csnk1a1
  suppression. Gsk3b row is Strong; others Moderate. S: Hyp-3B,
  T-Synergy.
- **Axis 3 / LR layer** (~8–10 rows): Gpc3_Unc5c; Entpd1_Adora1;
  L1cam_Ephb2; Gas6_Axl; Nrxn1_Nlgn1; Efna5_Epha4
  (broader, MG → Neuronal); L1cam_Cd9 (broader, Neuronal → MG);
  Sema6a_Plxna4 (sole 5/5-contrast G3 leader, sign-flips
  between amyloid-only and tau-on-amyloid). Mostly Moderate;
  Sema6a_Plxna4 is Strong. S: Hyp-3B, T-Synergy.
- **Axis 3 / pathway layer** (~4–6 rows): OXPHOS / ETC family (8
  GO terms); ribosomal compartment (5 GO CC + 1 GO MF at
  breadth = 4); Sun/Victor 2023 MG3_ribosome_biogenesis; DAM_up
  at interaction; MG-M3 (leader_score 23.4); all
  `interaction:3/mixed` cross-modality rows. Mostly Strong. S:
  Hyp-3B, T-Synergy.
- **Cross-axis Hyp-0 claims** (~3 rows): Cdk5 KSN footprint scaling;
  Rela sign-reversal across contrasts; Sema6a_Plxna4 single-LR
  divergence. Strong. S: Hyp-0 (+ secondary support to T-Tau-
  attenuates / T-Synergy as relevant).
- **Layer-level cross-axis observations** (~2–4 rows): TF
  asymmetric reading (strong / non-finding / strong); kinase
  asymmetric reading (same shape); LR refined-positive at axis 2;
  the layers AGREE that the synaptic axis is not regulator-driven
  at this resolution. Mostly Moderate. S: Hyp-0 (cross-axis pattern
  claim) + T-Compartment-suppression (the synaptic axis
  characterisation).

Expected confidence-grade distribution (target H2 outcome): Strong
≈ 12–18 rows; Moderate ≈ 18–30 rows; Suggestive ≈ 10–18 rows.
Expected support / contradict pattern (rough preview of H3): Hyp-1B
will outscore Hyp-1A (Rela sign-reversal + asymmetric reading);
Hyp-2B will outscore Hyp-2A on axis-restricted scope but the
cross-tool Pros1_Mertk row will narrow the margin; Hyp-3B will
outscore Hyp-3A strongly; T-Inflammation will receive substantial
support from axis-1 pathway rows but will be contradicted by
multiple axis-3 rows (the synergy contradicts pure-additivity);
T-Synergy will be the highest-support theme overall; Hyp-0 will
have 3 Strong rows.

#### H1 deliverable 3: H2 ledger source-priority order (locked here for execution)

H2 walks the source list in this order, emitting one claim row per
atomic finding. The order is most-distilled-first to keep claims
anchored in named-verdict prose where possible:

1. Verdict TSVs (3 rows each, per-axis): `tf_activity_verdict.tsv`,
   `kinase_activity_verdict.tsv`, `ccc_lr_verdict.tsv`. Each row's
   `evidence_summary` paragraph names the top claims for the axis
   at the layer; H2 should emit one ledger row per atomic
   evidence_summary clause.
2. Axis-restricted TSVs (per-axis × per-source detail):
   `tf_activity_axis_restricted.tsv`,
   `kinase_activity_axis_restricted.tsv`,
   `ccc_lr_axis_restricted.tsv`. Supplement the verdict claims
   with axis × tool detail.
3. Cross-tool / cross-modality leader boards:
   `pathway_survey_unified_leaderboard.tsv` (144 rows; pathway
   claims),
   `tf_activity_unified_leaderboard.tsv` (TF cross-source rows),
   `ccc_lr_cross_tool_leaderboard.tsv` (LR cross-tool rows; the
   Pros1_Mertk nuance row source).
4. Integration / concordance:
   `integration_logfc_per_contrast.tsv`, `modality_concordance.tsv`.
5. Pathway / module support: `hdwgcna_module_de.tsv`,
   `hdwgcna_module_enrichment.tsv`, `hdwgcna_module_hubs.tsv`,
   `nebula_per_state_localisation.tsv`,
   `triangulation_summary.tsv`,
   `triangulation_top_interaction.tsv`.
6. Multi-tool LR detail: `ccc_multinichenet_*.tsv`,
   `ccc_cellchat_*.tsv`, `ifn_nlgf_asymmetry.tsv`.
7. Completed-plan outcome summaries (the structured prior reading):
   `storage/notes/completed/mechanism_layer_plan_2026-05-24.md`
   §"Outcome summary"; `storage/notes/completed/pathway_overhaul_plan_2026-05-23.md`
   §"Outcome summary". Use these to surface claims that were
   surfaced narratively in the prior session but might be missed by
   pure TSV walking.

H2 should commit a checkpoint TSV at row 40 to allow incremental
review before pushing on to the 80-row ceiling; H2 may decide to
stop at 50 if claim quality drops or to push beyond 80 if the
verdict TSVs surface additional atomic claims worth recording.

---

### Session H1 (original step body retained below for audit trail)

**Decision gate at start.** Two binary locks needed from the user via
`AskUserQuestion`: (a) the competing-model set; (b) the claims-ledger
schema + confidence-grade rules. Both are foundational to H2-H5; both
must be locked before H2 begins.

**(a) Competing-model set — default proposal:** three per-axis
contests with two models each (six model definitions), plus an open
cross-axis claim slot. The contests are:

- **Contest 1 (axis: amyloid_activation):**
  - **Hyp-1A — Tau-independent amyloid program.** Amyloid drives a
    canonical DAM-state inflammatory program; the program is
    qualitatively identical in MAPTKI and P301S backgrounds; tau
    contributes no inflammatory modulator. The interaction contrast
    at axis 1 is zero or near-zero.
  - **Hyp-1B — Tau attenuates amyloid-driven NF-κB.** Amyloid drives
    NF-κB-led inflammation in MAPTKI; tau on the amyloid background
    SUPPRESSES this program (Rela sign-reversal: positive at
    `nlgf_in_maptki` / `nlgf_in_p301s`, negative at `interaction`).
    The interaction contrast at axis 1 is biologically informative.

- **Contest 2 (axis: synaptic_suppression):**
  - **Hyp-2A — Classical complement-pruning hypothesis.** NLGF-driven
    synaptic suppression is mediated by canonical C1q / C3 / Mertk-
    Pros1 / Tyro3-Axl complement-tagging + receptor-mediated
    engulfment. Predicts these LR pairs at the top of axis-2 cell-cell
    communication rankings.
  - **Hyp-2B — TREM2 / APP-fragment-mediated clearance.** NLGF-driven
    synaptic suppression is mediated by TREM2-mediated clearance +
    APP-fragment uptake + synaptic-adhesion modulation (Apoe_Trem2,
    App_Cd74, Cd200_Cd200r1, Vcan_Cd44, Ncam1 ligands). Predicts
    these pairs at the top instead of classical complement.

- **Contest 3 (axis: interaction_metabolic):**
  - **Hyp-3A — No interaction-specific mechanism.** The interaction
    contrast is arithmetically the difference of single-insult arms;
    no qualitatively new mechanism emerges at the tau × amyloid
    interface. Mechanism-layer signals at the interaction are
    smaller-than-amyloid-alone in magnitude or absent.
  - **Hyp-3B — Distinct synergy mechanism (Gsk3b / Myc / axon-guidance).**
    At the tau × amyloid interface a qualitatively new mechanism
    emerges: Gsk3b activation (the strongest single-kinase signal of
    the project, padj 0.002, supported by Mapk14 / Cdk5 / Cdk1
    canonical tau-kinase activation), Myc / Creb1 / Tp53 / Jun
    coordinated biosynthetic suppression at the TF layer, and axon-
    guidance / synaptic-adhesion intercellular rewiring at the LR
    layer. The mechanism is NOT inferrable from either single-insult
    arm because all three layers carry interaction-contrast-specific
    sign or magnitude patterns.

- **Cross-axis claim slot:** Cdk5-as-integrator — Cdk5 is the only
  kinase in the top 5 at all three axes, and its KSN substrate
  footprint scales with the axis (81 sites at amyloid universe, 263
  at synaptic universe, 153 at interaction universe). Does this
  cross-axis Cdk5 signal constitute a project-wide integrator claim?
  This is a candidate **Hyp-0 — Cdk5 cross-axis integrator** that the
  ledger can adjudicate as a stand-alone claim, parallel to (not
  competing with) the three per-axis contests.

**Alternative (a):** a single flat 4-model set spanning all three
axes — e.g. T-Inflammation (amyloid + tau both inflammatory; additive
DAM amplification), T-Compartment-suppression (NLGF suppresses
synaptic compartment via complement), T-Tau-attenuates (Rela sign-
reversal at axis 1), T-Synergy (Gsk3b/Myc interaction-specific
mechanism at axis 3). This is closer to the user's suggested phrasing
in the original AskUserQuestion option D ("e.g. canonical complement
pruning vs TREM2/APP mediated; uniform suppression vs Gsk3b-led
synergy") but loses the clean per-axis contest structure. Default
prefers per-axis contests because each axis is a distinct biological
question, the verdict TSVs already organise around axes, and the per-
contest adjudication is cleaner to render.

**(b) Claims-ledger schema — default proposal:** 14 columns,
exhaustively documented above in the locked-decisions table. Saved as
`storage/results/biological_model_claims_ledger.tsv`. The schema is
designed so adjudication in H3 is a structured groupby on
`supports_models` / `contradicts_models` rather than free-text
parsing.

**Alternative (b):** a leaner 8-column schema (drop `corroborating_evidence`,
`n_replicates_in_layers`, `notes`, and split `effect_size` into
`numeric_value` + `metric_type` columns). The lean schema is faster
to populate but loses the audit trail needed to defend the
confidence-grade assignments at H4 render time when a reader hovers
over a Suggestive claim and asks "why suggestive?". Default prefers
the 14-column schema for the audit trail.

**Confidence-grade definitions — default proposal:**

- **Strong** = (a) effect reproduces across ≥ 2 modalities with
  consistent sign AND (b) reaches FDR < 0.10 in primary inference
  tool AND (c) has cross-axis or cross-layer corroboration named
  explicitly in the corroborating-evidence column. Examples: Gsk3b at
  interaction (padj 0.002, kinase layer, corroborates Myc TF
  suppression and axon-guidance LR rewiring at the same axis).
- **Moderate** = effect is significant in primary inference tool
  (FDR < 0.10) but lacks ≥ 2-modality replication OR lacks cross-
  layer corroboration. Examples: Apoe_Trem2 at axis-2 rank 9 in
  LIANA+ but without cross-tool ranking at the same axis-relevant
  contrasts (the G3-vs-G4 complementarity property).
- **Suggestive** = effect is in top-N of primary inference tool but
  does not reach significance OR has known structural caveats (e.g.
  Interpretation A entanglement at axes 1 + 2 for TF and kinase
  layers; single-cell n_cells_used = 1 LR cells). Examples: most TF
  / kinase top-5 at axis 2 (Interpretation A entanglement).

**Alternative (confidence-grades):** a binary Strong / Weak grading.
Simpler but loses the meaningful middle tier where most of the
project's findings sit. Default prefers the three-tier scheme.

**Action plan for H1 session (this session, in order):**

1. Read the verdict TSVs + pathway leader board + completed-plan
   outcome summaries to confirm the model defaults above are
   data-grounded (mostly done in the H1 model-framework drafting;
   document a structured sketch of claim categories in this
   completion note for handoff to H2).
2. Present the default model set + alternative (a) + default ledger
   schema + alternative (b) + default confidence-grade rules +
   alternative (single-tier) to the user via `AskUserQuestion` (one
   multi-question call covering all three locks).
3. Lock the chosen values into the Locked-decisions table; remove
   from the Open-questions table.
4. Write a brief sketch of the H2 claim categories (axis × layer
   matrix) for the next session's handoff.
5. chown / re-knit (no-op for the plan-only edits) / verify /
   commit / continue to H2 only if context budget permits.

**Inputs for H1.** All under `storage/results/` unless noted:

- Verdict TSVs: `tf_activity_verdict.tsv`,
  `kinase_activity_verdict.tsv`, `ccc_lr_verdict.tsv` (3 rows each,
  per-axis with `evidence_summary` paragraphs).
- Axis-restricted TSVs (per-axis x per-source detail):
  `tf_activity_axis_restricted.tsv`,
  `kinase_activity_axis_restricted.tsv`,
  `ccc_lr_axis_restricted.tsv`.
- Cross-modality / cross-tool leader boards:
  `pathway_survey_unified_leaderboard.tsv` (D2; 144 rows),
  `tf_activity_unified_leaderboard.tsv` (E3),
  `ccc_lr_cross_tool_leaderboard.tsv` (G3).
- Integration / concordance: `integration_logfc_per_contrast.tsv`,
  `modality_concordance.tsv`.
- Pathway / module support: `hdwgcna_module_de.tsv`,
  `hdwgcna_module_enrichment.tsv`, `hdwgcna_module_hubs.tsv`,
  `nebula_per_state_localisation.tsv`,
  `triangulation_summary.tsv`, `triangulation_top_interaction.tsv`.
- Multi-tool LR detail: `ccc_multinichenet_*.tsv`,
  `ccc_cellchat_*.tsv`, `ifn_nlgf_asymmetry.tsv`.
- Completed-plan outcome summaries (the structured prior reading):
  `storage/notes/completed/mechanism_layer_plan_2026-05-24.md`
  §"Outcome summary" (lines 3239-3464);
  `storage/notes/completed/pathway_overhaul_plan_2026-05-23.md`
  §"Outcome summary" (lines 1630-1810).

**Outputs from H1.** This session is plan-edit only (decision gate +
ledger-schema lock). No new R code, no new TSV, no Rmd modification.
The plan file is the sole deliverable, updated with the locked model
set + schema + confidence rules in the Locked-decisions table.

---

### Session H2: build the claims ledger [DONE 2026-05-25]

**Completion note (2026-05-25).** H2 emitted
`storage/results/biological_model_claims_ledger.tsv` (75 rows × 14
columns, 50,526 bytes) via `scripts/build_biological_model_ledger.R`
(63 KB; 758 lines of carefully indented `row()` constructor calls plus
helper, validation, and summary blocks). The script enforces five
schema invariants on every load — claim_id uniqueness, supports /
contradicts subset-of-ENTITY_IDS membership, confidence_grade ∈
{Strong, Moderate, Suggestive}, direction ∈ {+, −, mixed}, and locked
axis / layer enumerations — so the ledger cannot be silently corrupted
on future regenerations and the H3 adjudication arithmetic can groupby
on the model-membership columns without parsing artefacts.

The ledger's distribution lands close to the H1 sketch but slightly
heavier on Strong-graded rows than the H1 expected range (Strong=24
vs the predicted 12-18; Moderate=39 vs 18-30; Suggestive=12 vs
10-18). The over-shoot on Strong is by design rather than by drift:
during ledger drafting the prior session applied the three-criteria
Strong rule (≥2-modality consistent-sign + FDR<0.10 in primary tool +
explicit cross-axis or cross-layer corroboration) consistently, and
the project's cross-modality + cross-layer richness produces more
qualifying rows than the H1 estimate anticipated. The Suggestive
count lands inside its expected band, confirming that the
Interpretation-A-entangled and n_cells_used=1 rows are honestly
flagged rather than upgraded. Per-axis the row count is
amyloid_activation 28 / synaptic_suppression 17 /
interaction_metabolic 25 / cross_axis 5 — axis 2 is the lightest
because the TF+kinase non-finding rows (3+3) are necessarily fewer
than the axis-1 / axis-3 multi-tool rich sections, but the 9 LR rows
+ 2 pathway rows cover the Hyp-2A-vs-Hyp-2B contest comprehensively
(5 of the 6 TF/kinase non-finding rows are honest hypothesis-
aligned negative evidence supporting T-Compartment-suppression). Per-layer the row count is tf 16
/ kinase 14 / lr 28 / pathway 15 / cross_layer 2 — LR dominates
because the G3 cross-tool leader board + G4 axis-restricted view +
G2 LIANA+ + multinichenet/cellchat collateral all contribute, and
the H2-040 Pros1_Mertk PROTOCOL row (one row each for the axis-
restricted-contradiction view and the cross-tool-partial-support
view) is the locked-in nuance from H1's Open-questions resolution.

Per-entity support / contradict counts realise the H1 outcome
preview point-for-point. The per-axis contest verdicts (H3 will
compute these formally) read Hyp-1B 24 vs Hyp-1A 21 raw support but Hyp-1B
24 vs Hyp-1A 18 net support — Hyp-1A is contradicted by exactly 3 rows
(H2-007 per-contrast asymmetric reading; H2-028 ifn_nlgf_asymmetry
8× P301S-specific gene differential; H2-072 cross-axis Rela sign-
reversal). Hyp-2B 8 vs Hyp-2A −4 with the Hyp-2A contradicts surfacing from
6 axis-restricted LR rows (H2-035 / H2-036 / H2-037 / H2-038 /
H2-039 / H2-045) and the lone Moderate Hyp-2A support coming from the
H2-040 PROTOCOL cross-tool Pros1_Mertk row. Hyp-3B 27 vs Hyp-3A −26 is
almost unanimously decided: every axis-3 row that supports Hyp-3B
also contradicts Hyp-3A by mutual exclusion (the two models are
contradictories within the axis-3 contest by construction), and no
row supports Hyp-3A on the project's evidence. Cross-axis themes: T-Synergy 28 leads, then
T-Inflammation 26, T-Compartment-suppression 17, Hyp-0 7, T-Tau-
attenuates 5; T-Tau-attenuates is sparse because the sign-reversal
evidence is concentrated in a handful of cross-axis rows (H2-028
ifn_nlgf_asymmetry + H2-050 A0A979HLR9 axis-3 + H2-057 Sema6a_Plxna4
axis-3 + H2-072 Rela cross-axis + H2-073 Sema6a_Plxna4 cross-axis)
rather than spread across many per-axis rows. Hyp-0 Cdk5 cross-axis
integrator gets 7 supports concentrated in 5 Strong rows (H2-033
KSN-footprint expansion at axis 2 + H2-057 Sema6a_Plxna4 secondary
Hyp-0 anchor at axis 3 + H2-071 cross-axis Cdk5 presence at all three
axes + H2-072 cross-axis Rela sign-reversal covariant + H2-073
cross-axis Sema6a_Plxna4 sign-flip) plus 2 Moderate per-axis Cdk5
rows (H2-009 axis-1 + H2-054 axis-3), even stronger than the H1
prediction of "3 Strong rows" for Hyp-0.

**Plan-spec deviations and bugfixes (explicit).** Three deviations
from the H1 sketch / Open-questions defaults are worth recording.
(i) Axis-1 pathway / module rows came in at 6 rather than 7 (rows
H2-023..H2-028; the seventh DAM-context row sketched in H1 was
folded into H2-026 / H2-027 to avoid double-counting WAM / HAM /
AD2 against DAM_up). (ii) Axis-3 pathway / module rows came in at 7
rather than the H1 sketch's 5 (rows H2-064..H2-070; the additional
rows are the H2-068 triangulation_top_interaction Colec12 row and
the H2-069 / H2-070 MG-M2 and MG-M3 module-DE-at-tau_in_nlgf rows
that surfaced as Moderate / Strong evidence during drafting and
were worth atomising rather than collapsing). (iii) Cross-axis rows
landed at 5 rather than the H1 sketch's 3, adding H2-074
(TF + kinase agree on axis-2 non-finding) and H2-075
(TF + kinase + LR agree on axis-3 synergy) as cross-layer
observation rows because both meet the Strong threshold via
cross-layer (rather than cross-axis) agreement and serve as the
load-bearing cross-layer-corroboration anchors for the H3 / H4
arguments.

During the H2 validation pass this session, an audit subagent
catalogued every claim-ID cross-reference in the script's
`corroborating_evidence` and `notes` columns and surfaced 8
references that pointed at the wrong rows (a renumbering-drift
artefact from the cross_axis section being added after some of the
per-axis rows already cited "the Hyp-0 row" by a now-stale H2-049 /
H2-064 / H2-065 / H2-053 anchor). All 8 were patched in place:
H2-002 / H2-006 / H2-009 / H2-015 / H2-033 / H2-050 / H2-054 /
H2-057 now point at the correct cross-axis anchors (H2-071 Cdk5
Hyp-0 + H2-072 Rela sign-reversal + H2-073 Sema6a_Plxna4 sign-flip)
or the correct per-axis anchors (H2-015 now correctly cites the
axis-2 Apoe_Trem2 rows H2-035 / H2-036 instead of the unrelated TF
non-finding rows H2-029 / H2-030). A separate spot-check subagent
verified that 5 load-bearing numeric anchors (Gsk3b padj 0.00152;
Apoe_Trem2 axis-2 rank 9 / mean 0.961; Pros1_Mertk leader_score
15.4 × 6 cross-tool rows; Sema6a_Plxna4 sole leader_score 25.4 with
the 5/5-contrast summary string; MG-M3 leader_score 23.4 with
module-DE logFC −2.91 padj 0.074) reproduce EXACTLY from the
source TSVs — no drift in the load-bearing quantitative claims.
The script re-runs cleanly post-edit with all five validation
invariants intact.

**Files written.** `scripts/build_biological_model_ledger.R` (the
ledger builder, 758 lines) and
`storage/results/biological_model_claims_ledger.tsv` (the ledger
itself, 75 rows × 14 columns). Both chowned to rstudio:rstudio.
H2 wrote no Rmd modifications and no new helpers in R/ — the
ledger is consumed by H3 (adjudication) and H4 (chapter rendering)
sessions downstream, not in this session.

**Handoff to H3.** H3 walks the ledger row-by-row, computes per-
entity support / contradict / net_support arithmetic for all 11
entities (the values are previewed in this completion note above
and can serve as a cross-check on H3's emission), groups per-axis
contests into binary verdict rows (3 rows: amyloid_activation,
synaptic_suppression, interaction_metabolic), and emits both
`biological_model_adjudication.tsv` (11 entities × ~10 arithmetic
columns) and `biological_model_contest_verdicts.tsv` (3 contest
rows). H3 should build a small adjudication helper in
`scripts/build_biological_model_adjudication.R` (parallel to this
session's builder) and validate that the per-entity arithmetic
matches the preview-counts above modulo any per-axis breakdown
columns that this session did not compute. The H1 plan's
Open-questions table already specifies the tie-break rules
(net_support → n_strong_supports → n_moderate_supports) and the
hybrid scheme (per-axis contests get binary verdicts; Hyp-0 + 4 themes
get stand-alone counts with axis-column breakdown).

---

### Session H2 (original step body retained below for audit trail)

**Goal.** Emit `storage/results/biological_model_claims_ledger.tsv`
with one row per atomic mechanism claim, populated per the H1-locked
schema. Target row count: 40-80 rows. The ledger is the structured
evidence base for H3 adjudication and the centrepiece of the H4
chapter.

**Per-claim extraction protocol** (apply uniformly):

1. Identify the atomic claim (one sentence, one biological assertion,
   one direction).
2. Locate the primary evidence source (the most distilled TSV / row
   that names the effect; usually a verdict TSV or the pathway leader
   board).
3. Locate corroborating evidence (other TSVs / sections that name
   the same effect from a different angle); semicolon-join the
   references in the corroborating-evidence column.
4. Count replicates: how many modalities show consistent-sign
   significance? How many mechanism layers? Populate the two
   `n_replicates_*` columns.
5. Assign the confidence grade per the H1-locked rules.
6. Adjudicate per model: which of Hyp-1A/Hyp-1B/Hyp-2A/Hyp-2B/Hyp-3A/Hyp-3B/Hyp-0 does
   this claim support, and which does it contradict? Populate the
   two `supports_models` / `contradicts_models` columns. A single
   claim can support one model and contradict another at the same
   axis (e.g. an Apoe_Trem2 ranking supports Hyp-2B and contradicts Hyp-2A).
7. Populate the `notes` column with caveats (Interpretation A
   entanglement, n_cells_used = 1, sign-reversal, opposite-direction
   companion, etc.).

**Anticipated claim categories** (rough sketch from H1 inspection of
the verdict TSVs; H2 finalises):

- **Axis 1 / TF layer:** Spi1 activation; Nfkb1 activation; Sp3
  activation; Rreb1 repression; A0A979HLR9 (CollecTRI complex
  source) activation; Rela activation at NLGF contrasts;
  Rela suppression at interaction contrast (Hyp-1B signature).
- **Axis 1 / kinase layer:** Cdk2 activation; Cdk5 activation
  (p301s-biased); Mapk8 activation (p301s-biased); Csnk1e activation;
  Camk2g repression.
- **Axis 1 / LR layer:** Adam11_Itga4 / Gas6_Axl / L1cam_Egfr /
  Omg_Rtn4r / Rspo3_Lgr4 top-5 (axon-guidance bias); Cd200_Cd200r1
  (rank 6); Apoe_Trem2 (rank 9); App_Cd74 (rank 11); Tgfb1_Sdc2
  (rank 34); 55/100 top-100 microglia-involving.
- **Axis 1 / pathway layer:** MG-M2 module (cross-modality
  consistent-sign = 4 at nlgf_in_p301s); GOBP_ADAPTIVE_IMMUNE_RESPONSE;
  MHC II antigen-presentation family; custom DAM_up / WAM /
  HAM_disease / AD2.
- **Axis 2 / TF layer:** Honest non-finding (M-shared claim);
  Interpretation A entanglement caveat.
- **Axis 2 / kinase layer:** Honest non-finding (M-shared);
  Interpretation A entanglement; Cdk5 KSN footprint expansion
  (81 → 263 sites).
- **Axis 2 / LR layer:** Apoe_Trem2 (MG_DAM → MG_DAM rank 9);
  Apoe_Trem2 (Astrocyte → MG_DAM rank 10); App_Cd74 (Oligodendrocyte
  → Microglia rank 12); 49/100 top-100 microglia-involving;
  6,459 axis-2 cells (largest axis); C1qb_Lrp1 (rank 775, the
  first complement pair); no C3 / Mertk-Pros1 / Tyro3-Axl in top
  800; Hyp-2B supported, Hyp-2A contradicted.
- **Axis 2 / pathway layer:** Six synaptic GO CC terms at leader_score
  11.6-15.4 (consistent-sign negative at NLGF contrasts).
- **Axis 3 / TF layer:** Myc suppression; Creb1 suppression;
  Tp53 suppression; Jun suppression; A0A979HLR9 (interaction
  context, negative; sign-reversal vs axis 1).
- **Axis 3 / kinase layer:** Gsk3b activation (padj 0.002; the
  strongest single-kinase result; per-contrast amyloid alone −1.66,
  tau alone −1.47, interaction +3.79 to +4.32; **synergy IS the
  activation**); Mapk14 activation; Cdk5 activation (cross-axis);
  Cdk1 activation; Csnk1a1 suppression.
- **Axis 3 / LR layer:** Gpc3_Unc5c; Entpd1_Adora1; L1cam_Ephb2;
  Gas6_Axl; Nrxn1_Nlgn1; Efna5_Epha4 (broader, MG → Neuronal);
  L1cam_Cd9 (broader, Neuronal → MG); Sema6a_Plxna4 (sole
  5/5-contrast G3 leader, sign-flips between amyloid-only and
  tau-on-amyloid).
- **Axis 3 / pathway layer:** OXPHOS / ETC family (8 GO terms);
  ribosomal compartment (5 GO CC + 1 GO MF at breadth=4);
  Sun/Victor 2023 MG3_ribosome_biogenesis; DAM_up at interaction;
  MG-M3 (leader_score 23.4); all `interaction:3/mixed` cross-modality.
- **Cross-axis claims (Hyp-0 slot):** Cdk5 KSN footprint scaling;
  Rela sign-reversal across contrasts; Sema6a_Plxna4 single-LR
  divergence.
- **Layer-level cross-axis observations:** TF asymmetric reading
  (strong / non-finding / strong); kinase asymmetric reading (same
  shape); LR refined-positive at axis 2; the layers AGREE that the
  synaptic axis is not regulator-driven at this resolution.

**Inputs for H2.** Same as H1's Inputs list above; reading is
selective per claim.

**Outputs from H2.**

- `storage/results/biological_model_claims_ledger.tsv` (the ledger
  itself, ~40-80 rows × 14 columns).
- Plan-spec update: H2 status DONE; completion-note paragraphs
  documenting row count, claim-category distribution, confidence-
  grade distribution, support / contradict counts per model
  (preview of H3), and any plan-spec deviations.

---

### Session H3: build the model adjudication table [DONE 2026-05-25]

**Completion note (2026-05-25).** H3 emitted both H3-locked TSVs via
`scripts/build_biological_model_adjudication.R` (15 KB; ~250 lines of
base-R groupby arithmetic with the same defensive-validation pattern
as the H2 builder). The script loads
`storage/results/biological_model_claims_ledger.tsv`, re-validates the
five schema invariants on load (claim-ID uniqueness, axis enumeration,
confidence-grade enumeration, supports / contradicts subset-of-
ENTITY_IDS membership for each of the 75 rows), builds two
75-row × 11-entity logical membership matrices via
`mat_logical("supports_models")` and `mat_logical("contradicts_models")`,
and reduces them with `colSums`-by-slice for the per-entity, per-axis,
and per-grade aggregates. The reduction is a single-pass groupby that
needs the matrices only once and avoids any string re-parsing during
arithmetic; this is fast (sub-second on the 75-row ledger) and
reproducible across regenerations.

**Outputs land at the expected paths and sizes.**
`storage/results/biological_model_adjudication.tsv` is 11 rows ×
**19 columns** (the H3 plan spec sketched "~10 arithmetic columns"; the
realised count is 19 because the per-axis × per-entity wide breakdown
adds 8 columns — support_at + contradict_at for each of the 4 axes —
plus 4 catalogue columns at the front (entity_id, entity_type,
contest_id, entity_name). The total of 4 catalogue + 7 overall-
arithmetic + 8 axis-wide = 19 lines up with the plan body's locked
arithmetic spec, just with a more generous catalogue header).
`storage/results/biological_model_contest_verdicts.tsv` is 3 rows ×
**16 columns** (contest_id + model_a × 6 columns including
support / contradict / net / n_strong / id / name + model_b × 6 same +
favoured_model + favoured_by_margin + tie_break_rule_invoked). Both
files are tab-separated, unquoted, empty-cell-as-NA, and chowned to
rstudio:rstudio.

**Headline arithmetic reproduces the H2 preview point-for-point.** The
H2 completion note previewed Hyp-1B 24 vs Hyp-1A 21 raw support / Hyp-1B 24 vs
Hyp-1A 18 net support; the H3 emission realises Hyp-1A support=21
contradict=3 net=18, Hyp-1B support=24 contradict=0 net=24 — exact match.
Hyp-2B 8 vs Hyp-2A −4 (H2 preview): Hyp-2A support=2 contradict=6 net=-4, Hyp-2B
support=8 contradict=0 net=8 — exact match. Hyp-3B 27 vs Hyp-3A −26 (H2
preview): Hyp-3A support=0 contradict=26 net=-26, Hyp-3B support=27
contradict=0 net=27 — exact match. Cross-axis theme leaderboard order
(H2 preview): T-Synergy 28 > T-Inflammation 26 > T-Compartment-
suppression 17 > T-Tau-attenuates 5, Hyp-0 7 — exact match. The H2
preview's "Hyp-0 gets 7 supports concentrated in 5 Strong rows" decomposes
in the H3 emission as `n_strong_supports = 5`, `n_moderate_supports = 2`,
`n_suggestive_supports = 0` — exact match (and stronger than the H1
prediction of "3 Strong rows" for Hyp-0). The H2 preview noted Hyp-2A's lone
Moderate support coming from the H2-040 PROTOCOL cross-tool
Pros1_Mertk row — H3 realises Hyp-2A `n_strong_supports = 1`,
`n_moderate_supports = 1`, `n_suggestive_supports = 0`, matching the
expected single Moderate row (the single Strong support is H2-052
classical-complement axis-3 cross-axis nuance).

**Per-axis contest verdicts are unambiguous; no tie-break rule was
invoked.** All three per-axis contests resolve at the first
arithmetic level (net_support) with comfortable margins: Hyp-1B beats
Hyp-1A by 6 (margin = 24 - 18); Hyp-2B beats Hyp-2A by 12 (margin = 8 - (-4));
Hyp-3B beats Hyp-3A by 53 (margin = 27 - (-26)). The `tie_break_rule_invoked`
column reads `no_tie_break_needed` for all three rows; the
`n_strong_supports` and `n_moderate_supports` tie-break logic is
plumbed and tested by hand but does not actually fire on this ledger.
Hyp-3B's margin of 53 is the largest of the three and reflects how
totally axis-3 falls to the synergy mechanism — every axis-3 row that
supports Hyp-3B also contradicts Hyp-3A by mutual exclusion (the two models
are contradictories within the axis-3 contest by construction), so
Hyp-3A's contradict_count of 26 mirrors Hyp-3B's support_count of 27 with
the off-by-one accounted for by H2-067 (the "all interaction:3/mixed
cross-modality rows" pathway claim that supports Hyp-3B but does not
explicitly contradict Hyp-3A because the row's `contradicts_models` cell
was left empty rather than asserting the strong-negative claim against
Hyp-3A; this is a deliberate H2 caveat). Hyp-1A's 3 contradictions come
from rows H2-007 (per-contrast asymmetric reading), H2-028
(ifn_nlgf_asymmetry 8× P301S-specific gene differential), and H2-072
(cross-axis Rela sign-reversal) — exactly as previewed in the H2
completion note.

**Per-axis × per-entity wide breakdown surfaces a clean axis-anchoring
pattern in the cross-axis themes.** Each cross-axis theme draws the
bulk of its support from a single axis: T-Inflammation 24 of 26
supports at amyloid_activation (with 2 trickling in from
interaction_metabolic via the inflammatory-rebound rows); T-
Compartment-suppression 16 of 17 supports at synaptic_suppression
(with 1 from the cross-axis section); T-Synergy 25 of 28 supports at
interaction_metabolic (with 3 from the cross-axis section); T-Tau-
attenuates 1 / 0 / 2 / 2 split across axes (the most diffuse theme,
which is biologically what you'd expect since the sign-reversal
mechanism is named at axis 1 + 3 with cross-axis covariance). Hyp-0's
1 / 1 / 2 / 3 split is the most evenly distributed of any entity,
confirming the integrator framing — Cdk5 fires meaningfully at all
three axes plus the cross-axis section. **This is the data pattern
H4 17.3 will name as a "single-axis anchoring" observation in the
cross-axis theme paragraph**; the Hyp-0 evenness will be named as a
separate point. This per-axis × per-entity granularity is exactly the
information the wide-column emission was designed to surface, and
was the deciding factor for choosing the wide format over a long-
format groupby table during script design.

**Plan-spec deviations (explicit).** Three minor deviations from the
H3 plan-body sketch are worth recording. (i) The plan spec said
"~10 arithmetic columns"; the realised count is 19 because the per-
axis wide breakdown was retained as a wide format rather than long-
form. The wide format is locked by the H3 arithmetic spec
("Emitted as wide columns `support_at_amyloid_activation`,
`contradict_at_amyloid_activation`, etc."), so the deviation is in
counting, not in the spec — the realised count includes the wide
columns and the catalogue header columns that the "~10" sketch did
not explicitly enumerate. (ii) The verdict table has 16 rather than
~10 columns, expanded to carry both model_a + model_b details
side-by-side (id, name, support, contradict, net_support,
n_strong_supports) so the verdict reads self-contained at H4 render
time without joining back to the per-entity table. (iii) The script
also emits a 4-way axis breakdown including `cross_axis` as the
fourth axis, even though the cross_axis "axis" is a structural label
in the ledger rather than a biological axis; this is necessary to
account for the 5 cross-axis ledger rows (H2-071..H2-075) and the Hyp-0
cross-axis attributions. None of these deviations affect downstream
arithmetic; H4 will consume both TSVs verbatim.

**Files written.** `scripts/build_biological_model_adjudication.R`
(the adjudication builder, ~250 lines),
`storage/results/biological_model_adjudication.tsv` (11 entities ×
19 columns, 1,615 B), and
`storage/results/biological_model_contest_verdicts.tsv` (3 contests ×
16 columns, 765 B). All three chowned to rstudio:rstudio. H3 wrote
no Rmd modifications and no new helpers in R/ — the adjudication
TSVs are consumed by H4 (chapter rendering) downstream, not in this
session, so the knit is unchanged from H2's state.

**Handoff to H4.** H4 builds `rmd/16_biological_model.Rmd` as section
17 with five subsections (17.1 model definitions, 17.2 ledger as
DT::datatable, 17.3 adjudication as two sub-tables, 17.4 open
questions / testable predictions, 17.5 integrated model statement);
H4 also wires the new child reference into `analysis.Rmd` between
lines 132-135. H4 should consume:
- the H1 ENTITY_CATALOGUE for 17.1 (the entity_type / contest_id /
  entity_name fields are already in the adjudication TSV at columns
  1-4, so H4 can derive the 17.1 table from a single load of
  `biological_model_adjudication.tsv` plus the
  `key_prediction` field hand-authored from the plan body's H1
  deliverable-1 prose);
- `biological_model_claims_ledger.tsv` for 17.2 (raw load into
  `DT::datatable` with the 14 ledger columns; H4 should consider
  hiding `corroborating_evidence` and `notes` by default and
  surfacing them via a column-toggle for readability);
- `biological_model_adjudication.tsv` + `biological_model_contest_verdicts.tsv`
  for 17.3 (per-axis-contest verdict table first as 17.3.1; per-
  entity stand-alone table second as 17.3.2; the latter sorted by
  net_support descending — the script's `theme_block <- theme_block[order(-net_support, -n_strong_supports), ]`
  already does this sort so H4 just reads the TSV in row order);
- the H2 ledger's `notes` column for 17.4 (each Moderate / Suggestive
  row's notes name what would upgrade it; H4 distils these into a
  numbered list of testable predictions).
H4 should also note in 17.3 the **single-axis anchoring** pattern
(each theme draws bulk support from one axis) + the **Hyp-0 evenness**
(1/1/2/3 split across axes) as the two "notable surprises in the
per-axis × entity breakdown" the plan body's 17.3.2 paragraph calls
for.

---

### Session H3 (original step body retained below for audit trail)

**Goal.** Aggregate the H2 ledger by model and emit
`storage/results/biological_model_adjudication.tsv` with one row per
model, per-axis support / contradict counts, and a favoured-model
verdict per contest.

**Adjudication arithmetic** (locked, accounting for the hybrid scheme):

- Per-(entity) aggregate (applies to all 11 entities — six per-axis
  models, Hyp-0, four themes): `support_count`, `contradict_count`,
  `net_support = support_count − contradict_count`,
  `support_to_contradict_ratio`, `n_strong_supports`,
  `n_moderate_supports`, `n_suggestive_supports`.
- Per-(entity × axis) aggregate: same arithmetic restricted to claims
  whose `axis` matches. Emitted as wide columns
  `support_at_amyloid_activation`, `contradict_at_amyloid_activation`,
  etc. For cross-axis themes this surfaces the per-axis evidence
  distribution that a stand-alone count would hide.
- Per-contest favoured-model verdict (applies to the three per-axis
  contests only): the model with the higher net_support wins; ties
  broken by `n_strong_supports`; further ties broken by
  `n_moderate_supports`. Emitted as a separate `contest_verdict`
  table with one row per contest (3 rows).
- Hyp-0 + cross-axis themes (5 entities): adjudicated as stand-alone
  evidence counts; no head-to-head winner. Emitted in the main
  per-entity table with the same arithmetic columns; the
  `contest_id` cell is `NA` for these rows.

**Outputs from H3.**

- `storage/results/biological_model_adjudication.tsv` (per-entity
  table, 11 rows × ~10 arithmetic columns).
- `storage/results/biological_model_contest_verdicts.tsv` (per-
  contest favoured-model verdict table, 3 rows: amyloid_activation,
  synaptic_suppression, interaction_metabolic).
- Plan-spec update: H3 status DONE.

---

### Session H4: render `rmd/16_biological_model.Rmd` late chapter [DONE 2026-05-25]

**Completion note (2026-05-25).** H4 emitted
`rmd/16_biological_model.Rmd` (32 KB; ~470 lines of prose + 4 kable
renders + 1 `DT::datatable` render + 1 catalogue-build chunk) and
wired it into `analysis.Rmd` as a new child reference between the
existing `child-ccc-mechanism` (rmd/15) and `child-session` (rmd/99)
chunks. The new chapter renders as **section 17 "Integrated
biological model"** with subsections 17.1, 17.2, 17.3 (containing
17.3.1 and 17.3.2), 17.4, and 17.5, exactly per the plan-body
sketch; "Session info and results manifest" renumbers cleanly from
section 17 to section 18 under the parent's `number_sections: true`
without manual intervention. The knit completes with **zero
`class="error"` and zero `class="warning"`** in the rendered
`analysis.html` (29.9 MB on disk), so the chapter's introduction of
DT into the project's tech stack does not regress the existing 1-16
sections.

The chapter is the project's first use of `DT::datatable` (the prior
sections render every table through `knitr::kable` with
`results = "asis"`). The DT widget on subsection 17.2 carries the
75-row H2 ledger with `filter = "top"`, `pageLength = 15`, a Copy /
CSV export button row, and a Column-visibility toggle that hides
five audit-trail columns by default (`primary_evidence_source`,
`corroborating_evidence`, `n_replicates_in_modalities`,
`n_replicates_in_layers`, `notes`) and re-shows them on click. This
honours the plan-spec preference for surfacing the audit-trail
columns only on demand — a default-visible 14-column DT would have
been visually noisy. Live filtering on `confidence_grade`,
`supports_models`, axis, layer, or any other column works
out-of-the-box via DT's column search panes; no helper code is
required beyond the constructor call.

Each chapter subsection's quantitative anchors reproduce the H2
preview counts and H3 emission TSVs point-for-point. Subsection 17.1
renders all 11 entities in a single 6-column kable grouped by
structural_type (per-axis contest models → cross-axis integrator →
cross-axis themes); the entity_statement and key_prediction cells are
hand-authored paraphrases of the H1-deliverable-1 prose
(plan-file lines 161-269) so the chapter is self-contained without
forcing the reader to consult the plan file. Subsection 17.3.1
renders the per-axis contest verdicts (Hyp-1B beats Hyp-1A by margin 6;
Hyp-2B beats Hyp-2A by margin 12; Hyp-3B beats Hyp-3A by margin 53) plus three
per-contest summary paragraphs naming the load-bearing supporting
claims by ledger row ID and the contradicting rows where applicable;
no tie-break rule was invoked, all three contests resolved at the
first arithmetic level. Subsection 17.3.2 renders the per-entity
support tally (T-Synergy 28 leads; Hyp-3A -26 trails) plus a per-axis
support breakdown for the 5 cross-axis entities that surfaces the
**single-axis anchoring** pattern (T-Inflammation 92% at axis 1,
T-Compartment-suppression 94% at axis 2, T-Synergy 89% at axis 3) and
the **Hyp-0 evenness** (1/1/2/3 split, the most balanced distribution
of any entity) — exactly the two diagnostic patterns the H3
completion note flagged as 17.3.2's payload. Subsection 17.4 lists
five numbered testable predictions grouped by the structural slot
each addresses (axis-2 non-finding closure; Hyp-2A scope ambiguity;
Hyp-3A resurrection; T-Tau-attenuates lifting; Hyp-0 Cdk5-substrate
enrichment as the load-bearing test). Subsection 17.5 states the
integrated reading in three paragraphs — per-axis verdicts,
cross-axis integration, and a synthesis paragraph framing the
project's evidence as "tau attenuates amyloid NF-kB at axis 1; NLGF
suppresses the synaptic compartment via TREM2/APP at axis 2; a
qualitatively new Gsk3b/Myc/axon-guidance mechanism emerges at the
interaction; Cdk5 ties the three axes together as a candidate
cross-axis integrator". The synthesis is explicitly hedged as "the
project's currently best-supported integrated reading", not as an
authoritative scientific conclusion, in line with the H1 anti-
anchoring guardrails.

**Plan-spec deviations (explicit).** Three deviations from the H4
plan-body sketch are worth recording. (i) Subsection 17.1 uses a
single combined 11-row kable rather than the plan-body's "optionally
render the three per-axis contests as a separate small table" — the
combined table groups contest models adjacent in the row order
(Hyp-1A/Hyp-1B/Hyp-2A/Hyp-2B/Hyp-3A/Hyp-3B), which makes the head-to-head pairing
visually evident without needing a second table. (ii) Subsection
17.3.2 renders **two** tables (the entity tally + the per-axis
support breakdown) rather than the plan-body's single-table
suggestion; the second table is necessary to surface the per-axis
anchoring numbers (24/2/0/0 for T-Inflammation, etc.) that the
single-table view would not expose at glance. (iii) Subsection 17.4
lists 5 numbered predictions rather than walking every Moderate /
Suggestive row's `notes` cell individually — the grouped format
gives the reader a cleaner narrative entry point than 51 row-by-row
bullets would. None of these deviations affect the chapter's
structural payload; the plan-spec subsections all render at the
designated section numbers with the designated content.

**Pre-knit smoke test.** The chapter was smoke-tested via `Rscript -e`
before the knit (per the per-step session loop's instruction to test
new cache-reading code against live caches before the ~3-minute knit
runs). The smoke test verified: (a) all three TSVs load cleanly via
`readr::read_tsv`; (b) the hand-authored `entity_prose` tibble's 11
entity_ids match the 11 entity_ids in the adjudication TSV with no
missing or extra IDs on either side; (c) the entity-catalogue join
produces 11 rows in the expected structural-type order; (d) the
verdict and entity-tally views render with the expected row counts
and column shapes; (e) the per-axis anchor view reproduces the
single-axis anchoring pattern (T-Inflammation 24/0/0/0 at axis 1, etc.)
that the H3 completion note previewed; (f) the `DT::datatable`
constructor builds without errors. The knit's runtime added ~10
seconds to the total wall-clock (the DT widget is the only
non-instant render in the new chapter).

**Files written.** `rmd/16_biological_model.Rmd` (32 KB; the new
chapter) and one new child-reference line in `analysis.Rmd`. Both
chowned to rstudio:rstudio. `analysis.html` was re-knitted and
verified to render section 17 with zero errors and zero warnings.

**Handoff to H5.** H5 adds the upfront "Headline integrated model"
snippet to `analysis.Rmd` between the existing scientific-question
content and the data section's first child reference. The snippet is
2-3 short paragraphs (~150 words) sourced from subsection 17.5; the
plan spec is locked at 3 paragraphs (paragraph 1 names the three
per-axis contests + references section 17; paragraph 2 states the
favoured-model verdict per axis in one clause each + names the
cross-axis-theme leaders; paragraph 3 is optional and notes the
strongest single-mechanism result of the project (Gsk3b at padj
0.00152) and the headline refinement of a literature-canonical
hypothesis (TREM2/APP over classical complement at axis 2)). H5
should also move the plan file to
`storage/notes/completed/biological_model_plan_<YYYY-MM-DD>.md` with
a `## Outcome summary` section per the plan-spec close protocol.

---

### Session H4 (original step body retained below for audit trail)

**Goal.** Build the new section-17 chapter as 5 subsections + wire
it into `analysis.Rmd`.

**Chapter structure (renders as section 17 "Integrated biological
model"):**

- **17.1 Competing biological models being adjudicated.** Render the
  eleven entities from H1 in a single structured table grouped by
  structural type: `entity_id`, `entity_type` (contest_model / m0 /
  theme), `contest_id` (NA for Hyp-0 + themes), `axis_or_scope`,
  `entity_name`, `entity_statement`, `key_prediction`. ~3/4 page.
  Optionally render the three per-axis contests as a separate small
  table that surfaces the Hyp-1A-vs-Hyp-1B / Hyp-2A-vs-Hyp-2B / Hyp-3A-vs-Hyp-3B head-
  to-head structure explicitly.
- **17.2 Claims ledger.** Render the H2 ledger as a `DT::datatable`
  (filterable, sortable; users can click confidence-grade = "Strong"
  or any of the 11 entity IDs in `supports_models` /
  `contradicts_models` to slice; users can filter by axis or layer).
  The ledger is the centrepiece of the chapter. ~1-2 pages depending
  on row count.
- **17.3 Model adjudication.** Render the H3 adjudication outputs as
  TWO sub-tables:
  - **17.3.1 Per-axis contest verdicts** — the H3
    `biological_model_contest_verdicts.tsv` (3 rows; one per axis
    contest). Add a one-paragraph per-contest summary block naming
    the favoured model and the headline supporting claims (3
    paragraphs total).
  - **17.3.2 Cross-axis entity support tallies** — the H3
    `biological_model_adjudication.tsv` (11 rows; one per entity).
    Sort by net_support descending. Add a single short paragraph
    naming the highest-support theme + the Hyp-0 evidence count + any
    notable surprises in the per-axis × entity breakdown (e.g. a
    cross-axis theme that gets all its support from a single axis).
  ~1.5 pages total.
- **17.4 Open questions and testable predictions.** Numbered list:
  per claim graded Moderate or Suggestive in H2, what additional
  evidence would upgrade it to Strong (or downgrade it)? Per
  contradicted per-axis model, what additional evidence would
  resurrect it? Per low-support theme, what evidence pattern would
  raise its standing? ~1/2 page.
- **17.5 Integrated model statement.** A 2-3 paragraph synthesis
  combining the favoured-model verdict per contest into a single
  biological reading of the data, with explicit cross-axis-theme
  framing (e.g. "the integrated reading favours T-Tau-attenuates +
  T-Synergy + T-Compartment-suppression-as-non-complement-route over
  T-Inflammation"). This is the upfront-mirror source text. ~1/2
  page.

**Wiring.** Add a new child reference to `analysis.Rmd`:

```r
```{r child-biological-model, child = "rmd/16_biological_model.Rmd"}
```
```

inserted between the existing `child-ccc-mechanism` (rmd/15) and
`child-session` (rmd/99) references at lines 132-135 of analysis.Rmd.

**Outputs from H4.**

- `rmd/16_biological_model.Rmd` (new file).
- Updated `analysis.Rmd` (one new child reference).
- Re-knitted `analysis.html` (verify section 17 renders cleanly with
  zero `class="error"` and zero `class="warning"`).

---

### Session H5: add upfront mirror + close plan [DONE 2026-05-25]

**Completion note (2026-05-25).** H5 inserted a new `## Headline
integrated biological model` subsection into `analysis.Rmd` between
the existing "Section content lives in `rmd/`" paragraph and the
first `child-data` reference (immediately after the `show-params`
chunk). The snippet renders as **section 1.1** under section 1
"Scientific question and design", which is the natural top-level
subsection slot for a synthesis-mirror that the data section's
content does not need to see. Three paragraphs total, ~205 words
(slightly over the plan's "~150 words" target but within the
3-paragraph cap): paragraph 1 names the three per-axis contests + 4
themes + Hyp-0 slot and refers to section 17 for the structured ledger;
paragraph 2 states the favoured-model verdict per axis in one clause
each (Hyp-1B/Hyp-2B/Hyp-3B with margins 6/12/53) plus the cross-axis theme
ranking (T-Synergy 28 > T-Inflammation 26 > T-Compartment-suppression
17 > Hyp-0 7 > T-Tau-attenuates 5) and the Hyp-0-evenness note; paragraph
3 names the two headline findings (Gsk3b padj 0.00152 + TREM2/APP-
over-classical-complement refinement at axis 2). The knit completes
with **zero `class="error"` and zero `class="warning"`** in the
rendered `analysis.html`; section numbering renumbers cleanly
(section 1.1 added, sections 2-18 unchanged from H4).

**Plan-spec deviations (explicit).** One deviation. (i) The snippet
is ~205 words rather than the plan's "~150 words" target. The
over-shoot is by ~37%; the extra wordcount is concentrated in
paragraph 2's per-axis verdict clauses, which the plan body's
example used as one-clause-each but which read more clearly as
two-clause-each with the margin and the supporting-evidence quip
(e.g. "Hyp-2B beats Hyp-2A by margin 12, with a residual TAM-kinase
Pros1-Mertk route narrowly surviving at specific cross-tool nodes").
Trimming to 150 words would have either dropped the Pros1-Mertk
scope nuance (which the H1 Open-questions resolution flagged as a
load-bearing scope qualification) or dropped the cross-axis theme
ranking (which the plan body listed as paragraph 2 content). The
over-shoot is preferred to either omission. The snippet remains
short relative to the long-form section 17 chapter (~205 words vs
~3500 words), so the upfront-mirror character is preserved.

**Files written.** `analysis.Rmd` (the new `## Headline integrated
biological model` subsection plus its 3-paragraph snippet, inserted
before the `child-data` chunk). Plan file is moved to
`storage/notes/completed/biological_model_plan_2026-05-25.md` after
this commit per the plan-spec close protocol; the Outcome summary
below records the synthesis-arc deliverables for future-session
reference. `analysis.html` was re-knitted and verified to render
sections 1.1, 17 (with all 5 subsections), and 18 cleanly.

---

### Session H5 (original step body retained below for audit trail)

**Goal.** Add the "Headline integrated model" snippet to
`analysis.Rmd` near the top, re-knit, verify, move plan to
completed/.

**Upfront-mirror placement.** Insert as a new subsection between the
"Scientific question and design" content and the data section's
first child reference. The snippet is 2-3 paragraphs (~150 words),
inline in `analysis.Rmd` rather than as a new child Rmd (because
it is too short to justify a child file). Structure:

> ### Headline integrated biological model
>
> [Paragraph 1: name the three per-axis contests adjudicated by the
> project's evidence; refer the reader to section 17 for the
> structured ledger.]
>
> [Paragraph 2: state the favoured-model verdict per axis in one
> clause each; name the cross-axis integrators if any.]
>
> [Paragraph 3 (optional): note the strongest single-mechanism result
> of the project (Gsk3b at padj 0.002) and the headline refinement of
> a literature-canonical hypothesis (TREM2/APP over classical
> complement pruning at axis 2).]

**Outputs from H5.**

- Updated `analysis.Rmd` (the upfront-mirror subsection added).
- Re-knitted `analysis.html` (verify zero errors / warnings).
- Plan moved to
  `storage/notes/completed/biological_model_plan_<YYYY-MM-DD>.md`
  with a `## Outcome summary` section documenting per-axis
  favoured-model verdicts, ledger row count, confidence-grade
  distribution, and the next-plan pointer (likely "no successor;
  user-directed at next session" unless the user signals otherwise).

---

## Anti-anchoring guardrails (re-read every session)

These exist because LLMs (this agent included) drift toward the most
salient prior finding. Each session must enforce them:

- **Always** derive each claim from a named TSV row / verdict
  paragraph / completed-plan outcome summary; never from training-
  data bias toward canonical AD-microglia stories.
- **Always** present each contest with its two competing models
  side-by-side; never collapse a contest into a single model before
  the ledger and adjudication finish.
- **Always** preserve the axis-2 TF + kinase non-finding honestly in
  the ledger (claims that document the entanglement, the in-universe
  target-count differential, the LR refined-positive); never paper
  over the non-finding with axis-1 enthusiasm.
- **Always** grade confidence transparently using the H1-locked
  rules; never assign Strong without the explicit cross-modality +
  cross-layer corroboration check.
- **Always** state which claim contradicts which model when relevant
  (e.g. Apoe_Trem2 ranking contradicts Hyp-2A); never report only the
  positive support without the contradictory direction.
- **Always** acknowledge the Interpretation A entanglement at axes
  1 + 2 in the `notes` column of every TF / kinase claim that
  inherits from the axis-1 ranking.
- **Always** alphabetise hdWGCNA modules in any chapter renders;
  never place MG-M3 first (inherited from the mechanism-layer
  guardrails).
- **Never** reintroduce Hallmark gene sets at any step (inherited).
- **Never** add new computation in any H phase. The synthesis is
  "structure what is already there"; new caches / DE runs / tool
  calls are out of scope.
- **Never** present the integrated model statement as authoritative
  scientific conclusion. The framing is "the project's current
  best-supported integrated reading on the available evidence",
  with explicit hooks for what external replication or experimental
  follow-up would adjudicate.
- **Never** drop a model from the contest just because it is
  contradicted. A contradicted model is a result; the ledger
  documents it.

## How to mark this plan complete

When H5 is DONE, move this file to
`storage/notes/completed/biological_model_plan_<YYYY-MM-DD>.md`,
`chown rstudio:rstudio`, and add a final `## Outcome summary` section
recording: (a) the locked-models set; (b) the per-contest favoured-
model verdict; (c) the per-claim ledger row count + confidence-grade
distribution; (d) the integrated model statement reproduced verbatim;
(e) any cross-axis claims surfaced in the Hyp-0 slot; (f) the next-plan
pointer.

## Outcome summary

### (a) Locked-models set

The synthesis adjudicated **eleven biological entities** organised
into a hybrid scheme locked at H1 (storage/notes/biological_model_plan.md
§H1 deliverable 1 prior to the move):

- **Per-axis binary contests (6 models, 3 contests).** Contest 1 at
  axis amyloid_activation: **Hyp-1A** Tau-independent amyloid program vs
  **Hyp-1B** Tau attenuates amyloid-driven NF-kB. Contest 2 at axis
  synaptic_suppression: **Hyp-2A** Classical complement-pruning hypothesis
  vs **Hyp-2B** TREM2 / APP-fragment-mediated clearance. Contest 3 at axis
  interaction_metabolic: **Hyp-3A** No interaction-specific mechanism vs
  **Hyp-3B** Distinct synergy mechanism (Gsk3b / Myc / axon-guidance).
- **Cross-axis integrator (1 entity, no contest):** **Hyp-0** Cdk5
  cross-axis integrator.
- **Cross-axis themes (4 entities, no head-to-head contests):**
  **T-Inflammation** (additive DAM-amplification across axes;
  subsumes Hyp-1A + Hyp-2A), **T-Compartment-suppression** (synaptic
  compartment suppressed via classical complement; subsumes Hyp-2A),
  **T-Tau-attenuates** (tau modifies amyloid rather than amplifying
  it; subsumes Hyp-1B), **T-Synergy** (qualitatively new mechanism at
  the tau x amyloid interface; subsumes Hyp-3B).

### (b) Per-contest favoured-model verdicts

Per `storage/results/biological_model_contest_verdicts.tsv` (3 rows,
emitted by `scripts/build_biological_model_adjudication.R`):

| Contest | Favoured | Net-support margin | Tie-break invoked |
|---|---|---|---|
| amyloid_activation    | **Hyp-1B** | 6 (Hyp-1B 24 vs Hyp-1A 18 net)   | no_tie_break_needed |
| synaptic_suppression  | **Hyp-2B** | 12 (Hyp-2B 8 vs Hyp-2A -4 net)   | no_tie_break_needed |
| interaction_metabolic | **Hyp-3B** | 53 (Hyp-3B 27 vs Hyp-3A -26 net) | no_tie_break_needed |

All three contests resolve at the first arithmetic level
(net_support); none of the secondary tie-break rules
(n_strong_supports, n_moderate_supports) were invoked. Hyp-3B's
margin of 53 is the largest and reflects how decisively the axis-3
synergy mechanism dominates the alternative.

### (c) Per-claim ledger row count + confidence-grade distribution

The H2 ledger at `storage/results/biological_model_claims_ledger.tsv`
records **75 atomic claims × 14 columns** (50,526 bytes; emitted by
`scripts/build_biological_model_ledger.R`). Distribution:

- **Confidence grade.** Strong = 24 rows; Moderate = 39 rows;
  Suggestive = 12 rows. (Strong over-shoots the H1 expected band
  of 12-18 because the project's cross-modality + cross-layer
  richness produces more qualifying rows than H1 anticipated; the
  three-criteria Strong rule was applied consistently throughout.)
- **Axis.** amyloid_activation = 28 rows; synaptic_suppression = 17
  rows; interaction_metabolic = 25 rows; cross_axis = 5 rows.
- **Layer.** lr = 28 rows; tf = 16 rows; pathway = 15 rows; kinase
  = 14 rows; cross_layer = 2 rows.
- **Per-entity support / contradict counts** (also in
  `storage/results/biological_model_adjudication.tsv`):

| Entity                    | Type             | Sup | Con | Net | Strong / Mod / Sug |
|---|---|---|---|---|---|
| T-Synergy                 | theme            | 28  | 0   |  28 | 9 / 17 / 2  |
| Hyp-3B                       | contest_model    | 27  | 0   |  27 | 8 / 17 / 2  |
| T-Inflammation            | theme            | 26  | 0   |  26 | 6 / 16 / 4  |
| Hyp-1B                       | contest_model    | 24  | 0   |  24 | 7 / 15 / 2  |
| Hyp-1A                       | contest_model    | 21  | 3   |  18 | 6 / 14 / 1  |
| T-Compartment-suppression | theme            | 17  | 0   |  17 | 8 / 5 / 4   |
| Hyp-2B                       | contest_model    |  8  | 0   |   8 | 5 / 3 / 0   |
| Hyp-0                        | m0 (integrator)  |  7  | 0   |   7 | 5 / 2 / 0   |
| T-Tau-attenuates          | theme            |  5  | 0   |   5 | 3 / 1 / 1   |
| Hyp-2A                       | contest_model    |  2  | 6   |  -4 | 1 / 1 / 0   |
| Hyp-3A                       | contest_model    |  0  | 26  | -26 | 0 / 0 / 0   |

### (d) Integrated model statement (reproduced verbatim from section 17.5)

> The project's evidence currently favours the following integrated
> biological reading of how amyloid reshapes microglia under different
> tau backgrounds. Each clause is grounded in the per-axis contest
> verdict (17.3.1) and cross-axis theme tally (17.3.2); none of it is
> asserted beyond what the H2 ledger and H3 adjudication support.
>
> **Per-axis verdicts.** At the amyloid-activation axis (axis 1), the
> project favours **Hyp-1B over Hyp-1A** by net-support margin 6 — amyloid
> drives a robust NF-kB-led inflammatory program in MAPTKI, and tau
> on the amyloid background SUPPRESSES (rather than amplifies) that
> program at the interaction contrast, with Rela sign-reversal as
> the load-bearing signature. At the synaptic-suppression axis
> (axis 2), the project favours **Hyp-2B over Hyp-2A** on axis-restricted
> scope by net-support margin 12 — NLGF-driven synaptic suppression
> is mediated by TREM2-mediated clearance plus APP-fragment uptake
> plus synaptic-adhesion modulation, NOT by canonical C1q / C3
> complement-tagging, although the cross-tool Pros1-Mertk signal
> narrowly keeps a TAM-kinase efferocytosis route alive at specific
> microglia-microglia nodes. At the interaction axis (axis 3), the
> project favours **Hyp-3B over Hyp-3A** by net-support margin 53 — a
> qualitatively new mechanism emerges at the tau x amyloid interface
> that is not present in either single-insult arm, anchored by the
> project's strongest single-kinase result (Gsk3b activation, padj
> 0.00152) and reinforced by coordinated TF biosynthetic suppression
> (Myc / Creb1 / Tp53 / Jun) and axon-guidance LR rewiring at the LR
> layer.
>
> **Cross-axis integration.** Layered on top of the three per-axis
> verdicts, the cross-axis themes rank **T-Synergy (28 supports, 9
> Strong) > T-Inflammation (26 supports, 6 Strong) > T-Compartment-
> suppression (17 supports, 8 Strong) > T-Tau-attenuates (5 supports,
> 3 Strong)** as the project's currently best-supported cross-axis
> readings. The per-axis-anchoring analysis (17.3.2) flags that all
> three of the highest-support themes draw 89-94% of their support
> from a single axis, so the cross-axis framing is mostly a
> re-statement of the per-axis verdicts in theme-overlay language —
> each theme is the cross-axis-friendly framing of its dominant
> axis's favoured-model verdict. The genuine cross-axis claim is
> **Hyp-0 (Cdk5 cross-axis integrator)**, which is the most evenly
> distributed entity in the adjudication (1 / 1 / 2 / 3 split across
> axes and the cross-axis section, 5 Strong rows total). Hyp-0's
> strongest evidence is the Cdk5 top-5-at-all-three-axes ranking
> with substrate footprint that scales with the axis universe size;
> its testable prediction (Cdk5-substrate enrichment at every axis)
> is not yet directly tested by a claim row and is the project's
> highest-value follow-up.
>
> **Integrated reading.** Putting these together: in MAPTKI mice,
> amyloid drives a canonical DAM inflammatory program with NF-kB at
> its transcriptional core (Hyp-1A's matched-direction support is
> genuine; the program is real). On top of that, tau pathology in
> the P301S background reshapes how amyloid acts on microglia in
> two qualitatively distinct ways. First, tau ATTENUATES the
> amyloid-driven NF-kB inflammation specifically at the interaction
> (T-Tau-attenuates, anchored at axis 1's Hyp-1B sign-reversal
> verdict). Second, at the interaction itself, a qualitatively new
> mechanism emerges that is not the arithmetic sum of single-insult
> arms: Gsk3b activation, biosynthetic / translational TF
> suppression, axon-guidance LR rewiring, mitochondrial / ribosomal
> pathway reorganisation (T-Synergy, anchored at axis 3's Hyp-3B
> verdict). Independently, the synaptic compartment is suppressed
> by NLGF via a TREM2 / APP-fragment-mediated route rather than
> classical complement (T-Compartment-suppression, anchored at
> axis 2's Hyp-2B verdict, with a residual TAM-kinase route partially
> surviving at specific microglia-microglia nodes). Cdk5 ties the
> three axes together as a candidate cross-axis kinase integrator
> (Hyp-0), most plausibly because Cdk5 substrates participate in all
> three axis biologies (DAM inflammatory regulation, synaptic
> adhesion / engulfment, tau / cytoskeletal phosphorylation). This
> reading is what the project's evidence currently supports; it is
> NOT an authoritative scientific conclusion, and the 17.4 testable
> predictions name the experiments that would either harden, refine,
> or overturn it.

### (e) Cross-axis claims in the Hyp-0 slot

Hyp-0 (Cdk5 cross-axis integrator) accumulated **7 supports / 0
contradicts / net 7 / 5 Strong + 2 Moderate** across the ledger.
The five Strong-graded Hyp-0-supporting rows are:

- **H2-033** (kinase, axis 2) — Cdk5 KSN footprint expansion from 81
  amyloid-universe sites to 263 synaptic-universe sites (3.2x
  increase, reflecting Cdk5's documented role in synaptic phospho-
  regulation).
- **H2-057** (lr, axis 3) — Sema6a-Plxna4 secondary Hyp-0 anchor at
  axis 3 (sole 5/5-contrast G3 leader with sign-flip between
  amyloid-only and tau-on-amyloid).
- **H2-071** (cross_axis, cross_layer) — Cdk5 in the kinase top 5 at
  all three axes (ranks 2 / 2 / 4); the load-bearing Hyp-0 claim.
- **H2-072** (cross_axis, cross_layer) — Cross-axis Rela sign-reversal
  covariant with Hyp-0 (NF-kB regulation interacts with Cdk5-mediated
  phosphorylation networks).
- **H2-073** (cross_axis, cross_layer) — Cross-axis Sema6a-Plxna4
  sign-flip, axis-3-anchored but with cross-layer corroboration via
  the kinase and pathway layers at the same axis.

The two Moderate Hyp-0 supports are per-axis Cdk5 rows: **H2-009** at
axis 1 (Cdk5 activation p301s-biased: maptki +0.50 vs p301s +3.00 at
the kinase layer) and **H2-054** at axis 3 (Cdk5 activation +2.10 at
the interaction). Together the 7 Hyp-0 rows cover all three axes plus
the cross-axis section, making Hyp-0 the most evenly distributed entity
in the adjudication (per-axis support split 1 / 1 / 2 / 3) and
vindicating the H1 decision to admit Hyp-0 as a stand-alone integrator
slot rather than forcing it into a per-axis contest.

### (f) Next-plan pointer

**No immediate successor plan.** Phase H closes the integrated-
biological-model synthesis arc; the project's three-axis cross-
modality verdicts (Phase D), three mechanism-layer verdicts (Phase
G), and integrated biological-model adjudication (Phase H) now read
together as a coherent end-to-end synthesis. The next session is
**user-directed at the next session start** — the canonical entry
point is a fresh-session survey-and-propose pass per the agent
launch protocol (storage/notes/reusable_prompt.md §"Step 2 (only if
no active plan exists)"). Candidate next directions that the H1-H5
arc has surfaced but not pursued, for the user's consideration:

- **Cdk5-substrate-restricted GSEA at all three axes.** Directly
  tests the Hyp-0 integrator framing (Cdk5-substrate enrichment
  expected at every axis if Hyp-0 is correct). The H1 → H4 sketch
  named this as the project's highest-value follow-up.
- **Pros1-Mertk cell-type-resolved validation in GeoMx microglia
  AOIs.** Settles the Hyp-2A scope ambiguity (whether the cross-tool
  Pros1-Mertk signal is genuinely microglia-microglia or driven by
  a different sender).
- **Per-substate NEBULA on MG_homeostatic at the interaction
  contrast.** Tests whether NF-kB attenuation is global or
  restricted to activated states (would lift T-Tau-attenuates out
  of its current 5-row position).
- **Permutation null on the axis-3 gene universe.** Tests whether
  the interaction-restricted gene universe is a poor proxy for true
  interaction biology (would weaken or strengthen Hyp-3B accordingly).
- **Manuscript draft.** The synthesis is now structurally complete;
  the document at `analysis.html` carries the full IMRAD-style
  narrative in 18 sections + section 1.1 upfront mirror. A
  manuscript draft (Methods + Results + Discussion + Figures) would
  be the natural next deliverable if the user is ready to convert
  the analysis into a publishable paper.

The five candidates above are non-exhaustive; the user may choose
any of them, propose a different direction, or sign-off on the
synthesis as currently complete. None of the candidates is
pre-locked — the next session begins with the survey-and-propose
protocol unless the user signals otherwise.
