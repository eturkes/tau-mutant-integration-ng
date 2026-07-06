#!/usr/bin/env python3
"""Inventory report prose blocks and planned visual replacements.

Counts human-facing Quarto text while skipping YAML, executable code bodies,
ordinary source comments, and generated output. Quarto captions/alt metadata
(`fig-cap`, `tbl-cap`, `fig-alt`) are retained as report-facing text.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from html.parser import HTMLParser
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


DEFAULT_QMDS = (
    "index.qmd",
    "_qc.qmd",
    "_story.qmd",
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
    ("_qc.qmd", "Quality control"): "collapsed-qc-audit",
    ("_qc.qmd", "Modalities loaded"): "qc-modality-table",
    ("_qc.qmd", "snRNAseq microglia: design balance"): "fig-qc-genotype-batch",
    ("_qc.qmd", "snRNAseq microglia: quality metrics"): "fig-qc-depth;fig-qc-fractions",
    ("_qc.qmd", "Sanity bounds (enforced)"): "collapsed-qc-audit",
    ("_story.qmd", "Scientific story"): "fig-story-core;fig-story-mechanism-crossmodality",
    ("_story.qmd", "Core evidence"): "fig-story-core",
    ("_story.qmd", "Mechanism and integration"): "fig-story-mechanism-crossmodality",
    ("_microglia.qmd", "Microglia: amyloid drives a homeostatic-to-DAM programme"): "fig-microglia-summary-board",
    ("_microglia.qmd", "Substate landscape"): "fig-microglia-substate-markers;fig-microglia-umap-substate;fig-microglia-score-triptych",
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
    ("_mechanism.qmd", "TF Activity"): "fig-mechanism-tf-interaction;fig-mechanism-tf-focus",
    ("_mechanism.qmd", "NF-kB Attenuation"): "fig-mechanism-nfkb-discordance",
    ("_mechanism.qmd", "Gsk3b And Kinase Activity"): "fig-mechanism-kinase-heatmap",
    ("_mechanism.qmd", "Mechanism status"): "fig-mechanism-status-board",
    ("_crossmodality.qmd", "Cross-modality: spatial, bulk, and divergence checks"): "fig-crossmodality-status-board",
    ("_crossmodality.qmd", "GeoMx Spatial DE"): "fig-crossmodality-geomx-volcano;fig-crossmodality-geomx-sensitivity",
    ("_crossmodality.qmd", "Bulk Proteome And Phospho"): "fig-crossmodality-bulk-run-index;fig-crossmodality-phospho-correction",
    ("_crossmodality.qmd", "Spatial Composition And Clearance Axis"): "fig-crossmodality-clearance-grid",
    ("_crossmodality.qmd", "Integrated Divergence"): "fig-crossmodality-symbol-matrix;fig-crossmodality-pathway-heatmap",
    ("_crossmodality.qmd", "Cross-modality status"): "fig-crossmodality-status-board",
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


@dataclass(frozen=True)
class HtmlBlock:
    line: int
    kind: str
    text: str

    @property
    def words(self) -> int:
        return len(WORD_RE.findall(self.text))


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


class ReportHtmlParser(HTMLParser):
    """Small Quarto-main-content parser for visible caption-only blockers."""

    IGNORED_TAGS = {"script", "style", "noscript", "template", "head"}
    DISPLAY_NON_TEXT_TAGS = {"canvas", "figure", "iframe", "img", "svg", "table"}
    VOID_TAGS = {
        "area", "base", "br", "col", "embed", "hr", "img", "input", "link",
        "meta", "param", "source", "track", "wbr",
    }

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.stack: list[dict[str, object]] = []
        self.blocks: list[HtmlBlock] = []
        self._capture_stack: list[dict[str, object]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr = {k: v or "" for k, v in attrs}
        classes = set(attr.get("class", "").split())
        parent = self.stack[-1] if self.stack else {}
        in_main = bool(parent.get("in_main")) or (
            tag == "main" and attr.get("id") == "quarto-document-content"
        )
        ignored = bool(parent.get("ignored")) or tag in self.IGNORED_TAGS
        in_caption = bool(parent.get("in_caption")) or tag in {"figcaption", "caption"}
        output_display_depth = int(parent.get("output_display_depth", 0))
        if "cell-output-display" in classes:
            output_display_depth += 1
        output_stdout_depth = int(parent.get("output_stdout_depth", 0))
        if "cell-output-stdout" in classes:
            output_stdout_depth += 1

        node = {
            "tag": tag,
            "classes": classes,
            "in_main": in_main,
            "ignored": ignored,
            "in_caption": in_caption,
            "line": self.getpos()[0],
            "output_display_depth": output_display_depth,
            "output_stdout_depth": output_stdout_depth,
        }
        if tag not in self.VOID_TAGS:
            self.stack.append(node)

        if not in_main or ignored:
            return
        if output_display_depth > 0 and (
            tag in self.DISPLAY_NON_TEXT_TAGS or {"quarto-figure", "quarto-float"} & classes
        ):
            for capture in self._capture_stack:
                if capture["kind"] == "html_cell_output_display":
                    capture["non_text_child"] = True
        if tag == "p" and not in_caption:
            self._start_capture("html_paragraph", node["line"])
        elif tag == "table" and not in_caption:
            self._start_capture("html_table", node["line"])
        elif "cell-output-display" in classes:
            self._start_capture("html_cell_output_display", node["line"])
        elif "cell-output-stdout" in classes:
            self._start_capture("html_stdout_provenance", node["line"])

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        n_stack = len(self.stack)
        self.handle_starttag(tag, attrs)
        if len(self.stack) > n_stack and self.stack[-1].get("tag") == tag:
            self.stack.pop()

    def handle_endtag(self, tag: str) -> None:
        self._close_captures(tag)
        if self.stack:
            self.stack.pop()

    def handle_data(self, data: str) -> None:
        if not data or not self.stack:
            return
        current = self.stack[-1]
        if not current.get("in_main") or current.get("ignored"):
            return
        if current["tag"] in {"script", "style"}:
            return
        for capture in self._capture_stack:
            capture["text"].append(data)

    def _start_capture(self, kind: str, line: object) -> None:
        self._capture_stack.append({"kind": kind, "line": int(line),
                                    "text": [], "non_text_child": False})

    def _close_captures(self, tag: str) -> None:
        if not self._capture_stack:
            return
        closing = {
            "html_paragraph": "p",
            "html_table": "table",
            "html_cell_output_display": "div",
            "html_stdout_provenance": "div",
        }
        remaining: list[dict[str, object]] = []
        for capture in self._capture_stack:
            if closing.get(str(capture["kind"])) != tag:
                remaining.append(capture)
                continue
            text = normalize_space(" ".join(str(x) for x in capture["text"]))
            if text and not capture.get("non_text_child") and not self._allowed_html_capture(str(capture["kind"]), text):
                self.blocks.append(HtmlBlock(int(capture["line"]), str(capture["kind"]), text))
        self._capture_stack = remaining

    @staticmethod
    def _allowed_html_capture(kind: str, text: str) -> bool:
        # Plot widgets/image outputs sometimes place fallback labels in display containers.
        # The S1 blocker is text-only output, so a display container with only Figure refs is allowed.
        return kind == "html_cell_output_display" and re.fullmatch(r"Figure\s+\d+(?:\.\d+)?", text)


def normalize_space(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def parse_html(path: Path) -> list[HtmlBlock]:
    parser = ReportHtmlParser()
    parser.feed(path.read_text(encoding="utf-8"))
    return parser.blocks


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


def source_kind_counts(blocks: list[Block]) -> Counter[str]:
    return Counter(block.kind for block in blocks)


def html_kind_counts(blocks: list[HtmlBlock]) -> Counter[str]:
    return Counter(block.kind for block in blocks)


def print_kind_summary(blocks: list[Block], html_blocks: list[HtmlBlock]) -> None:
    print("source_kind\tblocks\twords")
    source_counts = source_kind_counts(blocks)
    source_kinds = sorted(set(source_counts) | {"caption", "paragraph", "table_row", "list_item"})
    for kind in source_kinds:
        count = source_counts.get(kind, 0)
        words = sum(block.words for block in blocks if block.kind == kind)
        print(f"{kind}\t{count}\t{words}")
    print()
    print("html_kind\tblocks\twords")
    html_counts = html_kind_counts(html_blocks)
    html_kinds = (
        "html_paragraph",
        "html_table",
        "html_cell_output_display",
        "html_stdout_provenance",
    )
    for kind in html_kinds:
        count = html_counts.get(kind, 0)
        words = sum(block.words for block in html_blocks if block.kind == kind)
        print(f"{kind}\t{count}\t{words}")


def source_blockers(blocks: list[Block]) -> list[Block]:
    return [block for block in blocks if not (block.kind.startswith("h") or block.kind == "caption")]


def print_strict_blockers(blocks: list[Block], html_blocks: list[HtmlBlock]) -> None:
    if blocks:
        print()
        print("STRICT SOURCE BLOCKERS")
        print("qmd\tline\tkind\twords\tsection\ttext")
        for block in blocks:
            text = normalize_space(clean_text(block.text))
            print(f"{block.qmd}\t{block.line}\t{block.kind}\t{block.words}\t{block.section}\t{text}")
    if html_blocks:
        print()
        print("STRICT HTML BLOCKERS")
        print("html\tline\tkind\twords\ttext")
        for block in html_blocks:
            text = normalize_space(block.text)
            print(f"index.html\t{block.line}\t{block.kind}\t{block.words}\t{text}")


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
    parser.add_argument("--strict", action="store_true",
                        help="fail if source contains visible blocks other than headings/captions")
    parser.add_argument("--html", type=Path,
                        help="also inspect rendered report HTML for visible body prose/tables/provenance")
    args = parser.parse_args()

    files = args.files or [Path(p) for p in DEFAULT_QMDS]
    missing = [str(path) for path in files if not path.exists()]
    if missing:
        parser.error("missing qmd file(s): " + ", ".join(missing))

    blocks = read_blocks(files)
    html_blocks = parse_html(args.html) if args.html else []
    print_summary(blocks)
    if args.strict:
        print()
        print_kind_summary(blocks, html_blocks)
        src_blockers = source_blockers(blocks)
        print_strict_blockers(src_blockers, html_blocks)
        total = len(src_blockers) + len(html_blocks)
        if total:
            print(f"\nFAIL - caption-only strict gate found {total} blocker(s)")
            return 1
        print("\nPASS - caption-only strict gate")
    if args.manifest:
        write_manifest(blocks, args.manifest)
    elif not args.summary_only:
        print()
        write_manifest(blocks, None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
