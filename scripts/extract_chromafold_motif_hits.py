#!/usr/bin/env python3
"""Extract the exact ChromaFold CTCF motif calls needed by this benchmark.

ChromaFold distributes a 50-bp motif track rather than the underlying motif
intervals.  A deletion whose length is not divisible by 50 cannot be represented
correctly by shifting those bins.  This utility extracts the official
AnnotationHub AH104727 motif calls, reproduces the released reference track, and
saves the motif-level intervals needed to construct the derivative chromosome.

The optional junction scan is deliberately fail-closed.  It scores every motif
window crossing the new junction with the same three JASPAR 2022 matrices and
requires a large margin below the weakest score retained in AH104727.  It is a
negative junction-motif audit, not a replacement for genome-wide FIMO.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pickle
import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np


CHROM = "chr3"
DELETE_START = 4_976_067
DELETE_END = 4_976_790
DELETE_LEN = DELETE_END - DELETE_START
DISPLAY_START = 4_803_501
DISPLAY_END = 5_144_387
OUTPUT_BP = 10_000
INPUT_BP = 4_010_000
FINE_BP = 50
RESOURCE_NAME = "hg38.JASPAR2022_CORE_vertebrates_non_redundant_v2"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def chromosome_slice(root: SimpleNamespace, chrom: str) -> tuple[int, int]:
    names = np.asarray(root.seqnames.values.astype(str))
    lengths = np.asarray(root.seqnames.lengths, dtype=np.int64)
    match = np.flatnonzero(names == chrom)
    if match.size != 1:
        raise RuntimeError(f"expected one {chrom} run, found {match.size}")
    run = int(match[0])
    begin = int(lengths[:run].sum())
    return begin, begin + int(lengths[run])


def expand_rle_slice(rle: SimpleNamespace, begin: int, end: int) -> np.ndarray:
    values = np.asarray(rle.values.astype(str))
    lengths = np.asarray(rle.lengths, dtype=np.int64)
    run_ends = np.cumsum(lengths)
    run_starts = run_ends - lengths
    output = np.empty(end - begin, dtype="U1")
    for value, run_start, run_end in zip(values, run_starts, run_ends, strict=True):
        left = max(int(run_start), begin)
        right = min(int(run_end), end)
        if left < right:
            output[left - begin : right - begin] = value
    return output


def rasterize(
    starts_1based: np.ndarray,
    ends_1based: np.ndarray,
    strands: np.ndarray,
    scores: np.ndarray,
    first_bin: int,
    last_bin: int,
) -> np.ndarray:
    """Reproduce ChromaFold's R script, including reverse-strand precedence."""
    size = last_bin - first_bin + 1
    forward = np.zeros(size, dtype=np.float64)
    reverse = np.zeros(size, dtype=np.float64)
    for strand_code, target in ((1, forward), (-1, reverse)):
        selected = np.flatnonzero(strands == strand_code)
        for index in selected:
            start = int(starts_1based[index])
            end = int(ends_1based[index])
            score = float(scores[index])
            candidate_first = max(first_bin, (start - FINE_BP) // FINE_BP)
            candidate_last = min(last_bin, end // FINE_BP)
            for bin_index in range(candidate_first, candidate_last + 1):
                # The released R code uses inclusive IRanges [50b, 50b+50].
                overlap = max(
                    0,
                    min(end, bin_index * FINE_BP + FINE_BP)
                    - max(start, bin_index * FINE_BP)
                    + 1,
                )
                if overlap >= 10:
                    target[bin_index - first_bin] = score
    # The released pickle is log1p(raw score), with the reverse track replacing
    # the forward value when both strands occupy the same bin.
    return np.log1p(np.where(reverse > 0, reverse, forward))


def parse_meme(path: Path) -> dict[str, tuple[np.ndarray, int]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    motifs: dict[str, tuple[np.ndarray, int]] = {}
    index = 0
    while index < len(lines):
        if not lines[index].startswith("MOTIF "):
            index += 1
            continue
        name = lines[index].split(None, 1)[1]
        index += 1
        while index < len(lines) and not lines[index].startswith(
            "letter-probability matrix:"
        ):
            index += 1
        if index == len(lines):
            raise RuntimeError(f"missing probability matrix for {name}")
        width = int(re.search(r"w=\s*(\d+)", lines[index]).group(1))
        nsites = int(re.search(r"nsites=\s*(\d+)", lines[index]).group(1))
        matrix = np.asarray(
            [[float(value) for value in lines[index + row + 1].split()] for row in range(width)],
            dtype=np.float64,
        )
        motifs[name] = (matrix, nsites)
        index += width + 1
    return motifs


def pwm_score(sequence: str, matrix: np.ndarray, nsites: int) -> float:
    base_index = {"A": 0, "C": 1, "G": 2, "T": 3}
    pseudocount = 0.1
    adjusted = (matrix * nsites + pseudocount * 0.25) / (nsites + pseudocount)
    return float(
        sum(math.log2(adjusted[row, base_index[base]] / 0.25) for row, base in enumerate(sequence))
    )


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotationhub-rdata", required=True, type=Path)
    parser.add_argument("--reference-motif-pickle", required=True, type=Path)
    parser.add_argument("--common-inputs", required=True, type=Path)
    parser.add_argument("--meme-motifs", required=True, type=Path)
    parser.add_argument("--out-npz", required=True, type=Path)
    parser.add_argument("--out-audit", required=True, type=Path)
    args = parser.parse_args()

    try:
        import rdata
    except ImportError as exc:
        raise RuntimeError("install the pure-Python 'rdata' package for this extraction") from exc

    converted = rdata.conversion.convert(rdata.parser.parse_file(args.annotationhub_rdata))
    if RESOURCE_NAME not in converted:
        raise RuntimeError(f"{RESOURCE_NAME} is absent from {args.annotationhub_rdata}")
    root = converted[RESOURCE_NAME]
    begin, end = chromosome_slice(root, CHROM)
    starts = np.asarray(root.ranges.start[begin:end], dtype=np.int64)
    ends = starts + np.asarray(root.ranges.width[begin:end], dtype=np.int64) - 1
    strand_text = expand_rle_slice(root.strand, begin, end)
    strands = np.where(strand_text == "+", 1, np.where(strand_text == "-", -1, 0)).astype(np.int8)
    metadata = root.elementMetadata.listData
    scores = np.asarray(metadata["score"][begin:end], dtype=np.float64)
    names = np.asarray(metadata["name"][begin:end]).astype(str)

    first_center = (DISPLAY_START // OUTPUT_BP) * OUTPUT_BP
    last_center = int(np.ceil(DISPLAY_END / OUTPUT_BP)) * OUTPUT_BP
    required_start = first_center - 2_000_000 - 100
    required_end = last_center - 2_000_000 + INPUT_BP + DELETE_LEN + 100
    keep = (ends >= required_start) & (starts <= required_end)
    starts = starts[keep]
    ends = ends[keep]
    strands = strands[keep]
    scores = scores[keep]
    names = names[keep]

    first_bin = required_start // FINE_BP
    last_bin = required_end // FINE_BP
    reconstructed = rasterize(starts, ends, strands, scores, first_bin, last_bin)
    with args.reference_motif_pickle.open("rb") as handle:
        released_sparse = pickle.load(handle)[CHROM][0, first_bin : last_bin + 1]
    if hasattr(released_sparse, "toarray"):
        released_sparse = released_sparse.toarray()
    released = np.asarray(released_sparse).reshape(-1).astype(np.float64)
    difference = np.abs(reconstructed - released)
    if reconstructed.shape != released.shape or float(difference.max()) > 1e-12:
        raise RuntimeError(
            "motif-interval rasterization does not reproduce the released reference track: "
            f"max_abs={float(difference.max())}"
        )

    common = np.load(args.common_inputs, allow_pickle=False)
    derivative_sequence = bytes(common["deletion_sequence"]).decode("ascii")
    common_start = int(common["start"])
    junction = DELETE_START - common_start
    motif_matrices = parse_meme(args.meme_motifs)
    official_minimum_scores = {
        name: float(np.asarray(root.elementMetadata.listData["score"])[
            np.asarray(root.elementMetadata.listData["name"]).astype(str) == name
        ].min())
        for name in motif_matrices
    }
    junction_scores: dict[str, dict[str, object]] = {}
    for name, (matrix, nsites) in motif_matrices.items():
        width = matrix.shape[0]
        best = (-math.inf, -1, "")
        for offset in range(junction - width + 1, junction):
            sequence = derivative_sequence[offset : offset + width]
            for strand, oriented in (("+", sequence), ("-", reverse_complement(sequence))):
                value = pwm_score(oriented, matrix, nsites)
                if value > best[0]:
                    best = (value, offset, strand)
        junction_scores[name] = {
            "best_approximate_fimo_log_odds": best[0],
            "best_derivative_start_0based": common_start + best[1],
            "strand": best[2],
            "weakest_official_retained_score": official_minimum_scores[name],
            "margin_below_weakest_official_score": official_minimum_scores[name] - best[0],
        }
    if any(item["margin_below_weakest_official_score"] <= 5.0 for item in junction_scores.values()):
        raise RuntimeError("novel-junction motif is too close to the official FIMO retention boundary")

    args.out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out_npz,
        chromosome=np.asarray(CHROM),
        start_1based=starts,
        end_1based=ends,
        strand_code=strands,
        score=scores,
        motif_name=names,
        required_start_1based=np.int64(required_start),
        required_end_1based=np.int64(required_end),
        annotationhub_id=np.asarray("AH104727"),
        annotationhub_rdata_sha256=np.asarray(sha256(args.annotationhub_rdata)),
        reference_motif_pickle_sha256=np.asarray(sha256(args.reference_motif_pickle)),
    )
    audit = {
        "status": "PASS",
        "purpose": "motif-level derivative reconstruction for the 723-bp BHLHE40 deletion",
        "annotationhub_id": "AH104727",
        "annotationhub_rdata_sha256": sha256(args.annotationhub_rdata),
        "reference_motif_pickle_sha256": sha256(args.reference_motif_pickle),
        "output_npz_sha256": sha256(args.out_npz),
        "chromosome": CHROM,
        "extracted_hit_count": int(starts.size),
        "reference_reconstruction": {
            "bins_checked": int(reconstructed.size),
            "nonzero_bins": int(np.count_nonzero(released)),
            "max_absolute_difference": float(difference.max()),
            "mismatched_bins_at_1e-12": int(np.count_nonzero(difference > 1e-12)),
        },
        "novel_junction_scan": {
            "result": "NO_RETAINED_CTCF_MOTIF",
            "method": (
                "all windows crossing the derivative junction were scored against the three "
                "JASPAR 2022 CTCF PWMs; every best score was more than 5 score units below "
                "the weakest motif retained in the official AH104727 FIMO resource"
            ),
            "motifs": junction_scores,
        },
    }
    args.out_audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
