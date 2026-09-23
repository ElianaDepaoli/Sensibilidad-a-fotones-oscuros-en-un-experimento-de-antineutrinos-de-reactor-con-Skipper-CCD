#!/usr/bin/env python3
"""Verify the local paper corpus, manifests, provenance, and HTML links."""

from __future__ import annotations

import csv
import hashlib
import json
from html.parser import HTMLParser
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PAPERS = ROOT / "papers_reactor_dark_photon"
HTML_FILES = [
    ROOT / "report_reactor_dark_photon.html",
    ROOT / "rfh_reactor_dark_photon_en.html",
    ROOT / "rfh_reactor_dark_photon_es.html",
]


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    pdfs = sorted(PAPERS.glob("*.pdf"))
    assert len(pdfs) == 8, f"expected 8 PDFs, found {len(pdfs)}"

    with (PAPERS / "manifest.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == len(pdfs), "manifest/PDF count mismatch"
    for row, pdf in zip(rows, pdfs, strict=True):
        assert sha256(pdf) == row["sha256"], f"SHA-256 mismatch: {pdf.name}"

    for html_file in HTML_FILES:
        parser = LinkParser()
        parser.feed(html_file.read_text(encoding="utf-8"))
        for href in parser.links:
            if href.startswith(("http://", "https://", "#")):
                continue
            assert (html_file.parent / href).exists(), f"broken link in {html_file.name}: {href}"

    numbers = json.loads((ROOT / "provenance/numbers.json").read_text(encoding="utf-8"))
    required = {"value", "statement", "produced_by", "from_scratch", "from_library", "choices"}
    assert numbers and all(set(record) == required for record in numbers.values())
    assert (ROOT / "provenance/claims.yaml").read_text(encoding="utf-8").startswith("claims:\n")

    print(f"PASS: {len(pdfs)} PDFs match the manifest")
    print(f"PASS: {len(HTML_FILES)} HTML files parsed; all local links resolve")
    print(f"PASS: {len(numbers)} quantitative provenance records are structurally complete")


if __name__ == "__main__":
    main()
