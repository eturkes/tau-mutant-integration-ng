# Tau mutant integration

Reproducible R/Quarto integration of snRNA-seq, GeoMx, proteome, and
phosphoproteome evidence across the MAPTKI/P301S x amyloid design. The one
user-facing artefact is `report/tau-mutant-integration.html`.

## Run

Fresh Debian environment:

```sh
scripts/bootstrap/sysdeps.sh
scripts/bootstrap/rv.sh
scripts/bootstrap/quarto.sh
rv sync
scripts/check.sh
```

`scripts/check.sh` rebuilds and validates the report target. Raw inputs live
outside Git under `storage/data/`.

## Layout

| Path | Role |
|---|---|
| `R/core/` | constants, I/O, design, utilities, environment provenance |
| `R/analysis/` | modality-specific processing and inference |
| `R/report/` | compact figure data, plotting, and render helpers |
| `sections/` | report fragments included by `index.qmd` |
| `assets/` | report theme and bundled fonts |
| `scripts/bootstrap/` | pinned environment/bootstrap installers |
| `rv/` | R package-manager lock and activation machinery |
| `storage/` | ignored inputs, caches, QA output, and targets store |
| `report/` | ignored, generated standalone HTML |

Several leading-dot/underscore names are contracts, not parallel source trees:

- `.agent/` = tracked project memory, roadmap, live map, and decision history.
- `.claude/` = shared Claude Code settings and `/session-prompt` command.
- `.serena/` = project language-server configuration and read exclusions.
- `_targets.R` = targets pipeline entrypoint; `_targets.yaml` redirects its
  generated store to `storage/targets/`.
- `_quarto.yml` = Quarto project configuration.

Project conventions and live wiring are in `CLAUDE.md`, `.agent/memory.md`, and
`.agent/map.md`. Historical plans and paths remain historical records rather
than current layout documentation.
