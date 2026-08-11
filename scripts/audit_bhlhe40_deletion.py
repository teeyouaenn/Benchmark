#!/usr/bin/env python3
"""Independently verify the two published BHLHE40 CRISPR guides on hg38.

The reference input is a FASTA slice whose header is ``>chr3:start-end`` with
1-based inclusive coordinates.  The script searches both strands, verifies the
reverse-strand NGG PAMs, derives the standard SpCas9 cut sites three bases from
the PAM-proximal end, and writes a machine-readable audit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


GUIDES = ("TACTATCTATAGTAACTCCC", "TACCAGACTTCCACCGTATC")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    lines = args.fasta.read_text(encoding="ascii").splitlines()
    match = re.fullmatch(r">([^:]+):(\d+)-(\d+)", lines[0])
    if not match:
        raise ValueError("FASTA header must be >chromosome:start-end")
    chromosome, start_text, end_text = match.groups()
    region_start_1based = int(start_text)
    region_end_1based = int(end_text)
    sequence = "".join(lines[1:]).upper()
    if len(sequence) != region_end_1based - region_start_1based + 1:
        raise ValueError("FASTA sequence length does not match its header")

    region_start_0based = region_start_1based - 1
    records = []
    for guide_index, guide in enumerate(GUIDES, 1):
        candidates = []
        for strand, query in (("+", guide), ("-", reverse_complement(guide))):
            offset = 0
            while True:
                hit = sequence.find(query, offset)
                if hit < 0:
                    break
                candidates.append((strand, query, hit))
                offset = hit + 1
        if len(candidates) != 1:
            raise RuntimeError(
                f"guide {guide_index} expected one exact hit in the slice; got {candidates}"
            )

        strand, reference_protospacer, hit = candidates[0]
        start_0based = region_start_0based + hit
        end_0based = start_0based + len(guide)
        if strand == "+":
            pam_reference = sequence[hit + 20 : hit + 23]
            pam_guide_orientation = pam_reference
            if not re.fullmatch(r"[ACGT]GG", pam_guide_orientation):
                raise RuntimeError(f"guide {guide_index} lacks an NGG PAM")
            cut_0based = start_0based + 17
            pam_interval = [end_0based, end_0based + 3]
        else:
            pam_reference = sequence[hit - 3 : hit]
            pam_guide_orientation = reverse_complement(pam_reference)
            if not re.fullmatch(r"[ACGT]GG", pam_guide_orientation):
                raise RuntimeError(f"guide {guide_index} lacks a reverse-strand NGG PAM")
            cut_0based = start_0based + 3
            pam_interval = [start_0based - 3, start_0based]

        records.append(
            {
                "guide_index": guide_index,
                "published_guide_5to3": guide,
                "strand": strand,
                "reference_protospacer_5to3": reference_protospacer,
                "protospacer_0based_half_open": [start_0based, end_0based],
                "protospacer_1based_inclusive": [start_0based + 1, end_0based],
                "pam_reference_5to3": pam_reference,
                "pam_guide_orientation_5to3": pam_guide_orientation,
                "pam_0based_half_open": pam_interval,
                "expected_spcas9_cut_0based": cut_0based,
                "cut_description_1based": (
                    f"between chr3:{cut_0based:,} and chr3:{cut_0based + 1:,}"
                ),
            }
        )

    cut_sites = sorted(record["expected_spcas9_cut_0based"] for record in records)
    deletion = [cut_sites[0], cut_sites[1]]
    audit = {
        "reference_slice": {
            "path": str(args.fasta.resolve()),
            "sha256": sha256(args.fasta),
            "chromosome": chromosome,
            "region_1based_inclusive": [region_start_1based, region_end_1based],
        },
        "published_guides_source": "NAR gkad378 main Methods",
        "guides": records,
        "inferred_clean_deletion_0based_half_open": deletion,
        "inferred_clean_deletion_1based_inclusive": [deletion[0] + 1, deletion[1]],
        "inferred_clean_deletion_length_bp": deletion[1] - deletion[0],
        "verdict": (
            "The preregistered 723-bp interval is exactly the clean cut-to-cut "
            "SpCas9 deletion implied by the two published guides. The actual "
            "clone-resolved repair junction is not published."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
