#!/usr/bin/env python3
"""Inventory report prose blocks and planned visual replacements.

Counts human-facing Quarto text while skipping YAML, executable code bodies,
ordinary source comments, and generated output. Quarto captions (`fig-cap`,
`tbl-cap`, `fig-alt`) are retained because they render as report text.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


DEFAULT_QMDS = (
    "index.qmd",
    "_synthesis.qmd",
    "_qc.qmd",
    "_microglia.qmd",
    "_trajectory.qmd",
    "_mechanism.qmd",
    "_crossmodality.qmd",
)

WORD_RE = re.compile(r"[A-Za-z0-9]+(?:[-'][A-Za-z0-9]+)?")
FENCE_RE = re.compile(r"^\s*````*")
HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)$")
LIST_RE = re.compile(r"^\s*(?:[-*+]|\d+[.)])\s+(.*)$")
CHUNK_LABEL_RE = re.compile(r"^\s*#\|\s*label:\s*(.+?)\s*$")
CAPTION_RE = re.compile(r"^\s*#\|\s*(fig-cap|tbl-cap|fig-alt):\s*(.+?)\s*$")
ATTR_RE = re.compile(r"\s*\{#[^}]+\}\s*$")
INLINE_R_RE = re.compile(r"`r\s+[^`]+`")
CODE_SPAN_RE = re.compile(r"`([^`]+)`")
LINK_RE = re.compile(r"\[([^\]]+)\]\([^)]+\)")
REF_RE = re.compile(r"@[A-Za-z0-9_.:-]+")
HTML_RE = re.compile(r"<[^>]+>")


SECTION_SLOTS = {
    ("index.qmd", "Overview"): "fig-report-spine-schematic",
    ("_synthesis.qmd", "Synthesis: the compact answer"): "fig-synthesis-visual-abstract",
    ("_qc.qmd", "Quality control"): "collapsed-qc-audit",
    ("_qc.qmd", "Modalities loaded"): "qc-modality-table",
    ("_qc.qmd", "snRNAseq microglia: design balance"): "fig-qc-genotype-batch",
    ("_qc.qmd", "snRNAseq microglia: quality metrics"): "fig-qc-depth;fig-qc-fractions",
    ("_qc.qmd", "Sanity bounds (enforced)"): "collapsed-qc-audit",
    ("_microglia.qmd", "Microglia: amyloid drives a homeostatic-to-DAM programme"): "fig-microglia-summary-board",
    ("_microglia.qmd", "Substate landscape"): "fig-microglia-umap-substate;fig-microglia-score-triptych",
    ("_microglia.qmd", "Amyloid expands the DAM compartment"): "fig-microglia-composition-shift;fig-microglia-unit-composition",
    ("_microglia.qmd", "Composition is tested across the genotype contrasts"): "fig-microglia-composition-forest;fig-microglia-composition-concordance",
    ("_microglia.qmd", "The amyloid-to-DAM expression programme"): "fig-microglia-whole-volcano",
    ("_microglia.qmd", "The tau-by-amyloid interaction is under-powered, not absent"): "fig-microglia-whole-volcano",
    ("_microglia.qmd", "Per-substate differential expression"): "fig-microglia-substate-audit;fig-microglia-substate-volcano",
    ("_microglia.qmd", "Caveats and provenance"): "collapsed-microglia-audit",
    ("_trajectory.qmd", "The tau-amyloid synergy adds DAM cells rather than advancing them"): "fig-trajectory-logic-board",
    ("_trajectory.qmd", "Amyloid drives a large shift along the activation axis"): "fig-trajectory-pseudotime-shift;fig-trajectory-pt-density",
    ("_trajectory.qmd", "The interaction is composition, not progression"): "fig-trajectory-decomposition;fig-trajectory-kitagawa-forest",
    ("_trajectory.qmd", "A per-cell model corroborates the position shift"): "fig-trajectory-logic-board",
    ("_trajectory.qmd", "Reconciliation, robustness and concordance"): "fig-trajectory-concordance;fig-trajectory-audit",
    ("_trajectory.qmd", "Caveats and provenance"): "collapsed-trajectory-audit",
    ("_mechanism.qmd", "Mechanism: Myc signal, NF-kB check, and bulk kinase support"): "fig-mechanism-status-board",
    ("_mechanism.qmd", "Pathway Survey"): "fig-mechanism-project-pathway;fig-mechanism-go-dotplot",
    ("_mechanism.qmd", "TF Activity"): "fig-mechanism-tf-lollipop",
    ("_mechanism.qmd", "NF-kB Attenuation"): "fig-mechanism-nfkb-discordance",
    ("_mechanism.qmd", "Gsk3b And Kinase Activity"): "fig-mechanism-kinase-heatmap",
    ("_mechanism.qmd", "Synthesis"): "fig-mechanism-status-board",
    ("_crossmodality.qmd", "Cross-modality: spatial, bulk, and divergence checks"): "fig-crossmodality-status-board",
    ("_crossmodality.qmd", "GeoMx Spatial DE"): "fig-crossmodality-geomx-volcano;fig-crossmodality-geomx-sensitivity",
    ("_crossmodality.qmd", "Bulk Proteome And Phospho"): "fig-crossmodality-bulk-run-index;fig-crossmodality-phospho-correction",
    ("_crossmodality.qmd", "Spatial Composition And Clearance Axis"): "fig-crossmodality-clearance-grid",
    ("_crossmodality.qmd", "Integrated Divergence"): "fig-crossmodality-symbol-matrix;fig-crossmodality-pathway-heatmap",
    ("_crossmodality.qmd", "P4 Synthesis"): "fig-crossmodality-status-board",
}


@dataclass(frozen=True)
class Block:
    qmd: str
    block_id: str
    line: int
    kind: str
    section: str
    label: str
    text: str

    @property
    def words(self) -> int:
        return len(WORD_RE.findall(clean_text(self.text)))


def clean_text(text: str) -> str:
    out = INLINE_R_RE.sub(" ", text)
    out = LINK_RE.sub(r"\1", out)
    out = REF_RE.sub(" ", out)
    out = CODE_SPAN_RE.sub(r"\1", out)
    out = HTML_RE.sub(" ", out)
    out = re.sub(r"\{#[^}]+\}", " ", out)
    out = out.replace("&nbsp;", " ")
    out = out.replace("\\(", " ").replace("\\)", " ")
    return out


def strip_quotes(text: str) -> str:
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        return text[1:-1]
    return text


def slug(text: str) -> str:
    cleaned = clean_text(text).lower()
    words = WORD_RE.findall(cleaned)
    return "-".join(words[:8]) or "block"


def parse_qmd(path: Path) -> list[Block]:
    lines = path.read_text(encoding="utf-8").splitlines()
    blocks: list[Block] = []
    in_yaml = bool(lines and lines[0].strip() == "---")
    in_code = False
    chunk_label = ""
    section = ""
    pending: list[str] = []
    pending_line = 0
    seq = 1

    def add(kind: str, line: int, text: str, label: str = "") -> None:
        nonlocal seq
        text = ATTR_RE.sub("", text).strip()
        if not text:
            return
        block_id = f"{Path(path).stem.lstrip('_')}-{seq:03d}-{kind}-{slug(text)}"
        blocks.append(Block(path.name, block_id, line, kind, section, label, text))
        seq += 1

    def flush() -> None:
        nonlocal pending, pending_line
        if pending:
            add("paragraph", pending_line, " ".join(s.strip() for s in pending))
            pending = []
            pending_line = 0

    for idx, raw in enumerate(lines, start=1):
        line = raw.rstrip()
        if in_yaml:
            if idx > 1 and line.strip() == "---":
                in_yaml = False
            continue

        if FENCE_RE.match(line):
            if in_code:
                in_code = False
                chunk_label = ""
            else:
                flush()
                in_code = True
                chunk_label = ""
            continue

        if in_code:
            label_match = CHUNK_LABEL_RE.match(line)
            if label_match:
                chunk_label = strip_quotes(label_match.group(1))
                continue
            cap_match = CAPTION_RE.match(line)
            if cap_match:
                add("caption", idx, strip_quotes(cap_match.group(2)), chunk_label)
            continue

        stripped = line.strip()
        if not stripped:
            flush()
            continue
        if stripped.startswith("<!--"):
            flush()
            continue
        if stripped.startswith((":::", "{{<", "```")):
            flush()
            continue
        if stripped.startswith("|") or re.match(r"^:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)+$", stripped):
            flush()
            add("table_row", idx, stripped)
            continue

        heading_match = HEADING_RE.match(line)
        if heading_match:
            flush()
            section = ATTR_RE.sub("", heading_match.group(2)).strip()
            add(f"h{len(heading_match.group(1))}", idx, section)
            continue

        list_match = LIST_RE.match(line)
        if list_match:
            flush()
            add("list_item", idx, list_match.group(1))
            continue

        if not pending:
            pending_line = idx
        pending.append(line)

    flush()
    return blocks


def disposition(block: Block) -> tuple[str, str]:
    slot = SECTION_SLOTS.get((block.qmd, block.section), "chapter-board")
    text_l = clean_text(block.text).lower()
    if block.kind.startswith("h"):
        return "keep", "section-nav"
    if block.kind == "caption":
        return "caption", block.label or slot
    if block.kind == "table_row":
        return "figure", slot
    if "caveat" in text_l or "provenance" in text_l or "warning" in text_l:
        return "collapsed_audit", slot
    if "spatialdecon" in text_l or "blocked" in text_l or "not called" in text_l:
        return "figure", slot
    if "synthesis" in block.section.lower() or block.qmd in {"index.qmd", "_synthesis.qmd"}:
        return "schematic", slot
    if "sanity" in block.section.lower() or "quality control" in block.section.lower() or "sanity" in text_l:
        return "collapsed_audit", slot
    if "this section" in text_l or "this chapter" in text_l:
        return "schematic", slot
    if block.kind == "list_item":
        return "collapsed_audit", slot
    return "figure", slot


def read_blocks(files: list[Path]) -> list[Block]:
    blocks: list[Block] = []
    for path in files:
        blocks.extend(parse_qmd(path))
    return blocks


def print_summary(blocks: list[Block]) -> None:
    by_qmd: dict[str, list[Block]] = defaultdict(list)
    for block in blocks:
        by_qmd[block.qmd].append(block)

    print(
        "qmd\tblocks\twords\theading_blocks\tprose_blocks\t"
        "non_keep_prose_blocks\tnon_keep_prose_pct\tnon_keep_words"
    )
    total = Counter()
    for qmd in DEFAULT_QMDS:
        chapter_blocks = by_qmd.get(qmd, [])
        words = sum(b.words for b in chapter_blocks)
        heading_blocks = [b for b in chapter_blocks if b.kind.startswith("h")]
        prose_blocks = [b for b in chapter_blocks if not b.kind.startswith("h")]
        non_keep_prose = [b for b in prose_blocks if disposition(b)[0] != "keep"]
        non_keep_words = sum(b.words for b in non_keep_prose)
        pct = 100 * len(non_keep_prose) / len(prose_blocks) if prose_blocks else 100
        print(
            f"{qmd}\t{len(chapter_blocks)}\t{words}\t{len(heading_blocks)}\t"
            f"{len(prose_blocks)}\t{len(non_keep_prose)}\t{pct:.1f}\t{non_keep_words}"
        )
        total.update(
            blocks=len(chapter_blocks),
            words=words,
            heading_blocks=len(heading_blocks),
            prose_blocks=len(prose_blocks),
            non_keep_prose=len(non_keep_prose),
            non_keep_words=non_keep_words,
        )
    total_pct = 100 * total["non_keep_prose"] / total["prose_blocks"] if total["prose_blocks"] else 100
    print(
        f"TOTAL\t{total['blocks']}\t{total['words']}\t{total['heading_blocks']}\t"
        f"{total['prose_blocks']}\t{total['non_keep_prose']}\t{total_pct:.1f}\t"
        f"{total['non_keep_words']}"
    )


def write_manifest(blocks: list[Block], path: Path | None) -> None:
    fieldnames = (
        "qmd",
        "block_id",
        "line",
        "kind",
        "section",
        "words",
        "disposition",
        "target_slot",
        "label",
        "text",
    )
    out = path.open("w", encoding="utf-8", newline="") if path else sys.stdout
    try:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for block in blocks:
            disp, slot = disposition(block)
            row = {
                "qmd": block.qmd,
                "block_id": block.block_id,
                "line": block.line,
                "kind": block.kind,
                "section": block.section,
                "words": block.words,
                "disposition": disp,
                "target_slot": slot,
                "label": block.label,
                "text": clean_text(block.text).strip(),
            }
            writer.writerow({k: ("." if v == "" else v) for k, v in row.items()})
    finally:
        if path:
            out.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", type=Path, help="QMD files to inventory")
    parser.add_argument("--manifest", type=Path, help="write block-level TSV manifest")
    parser.add_argument("--summary-only", action="store_true", help="skip manifest stdout when no --manifest is given")
    args = parser.parse_args()

    files = args.files or [Path(p) for p in DEFAULT_QMDS]
    missing = [str(path) for path in files if not path.exists()]
    if missing:
        parser.error("missing qmd file(s): " + ", ".join(missing))

    blocks = read_blocks(files)
    print_summary(blocks)
    if args.manifest:
        write_manifest(blocks, args.manifest)
    elif not args.summary_only:
        print()
        write_manifest(blocks, None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
