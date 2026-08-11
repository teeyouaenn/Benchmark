#!/usr/bin/env python3
"""Prepare immutable BHLHE40 WT/deletion inputs and experimental reference.

The primary deletion arm changes DNA only and holds the control K562 DNase
track fixed.  A second, explicitly target-assisted surrogate maps the deletion
Hi-TrAC endpoint profile onto the empirical distribution of the control DNase
track.  The latter must never be reported as de novo prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pyBigWig
from pyfaidx import Fasta
from scipy.stats import rankdata, spearmanr


CHROM = "chr3"
START = 3_000_000
END = 7_000_000
DELETE_START = 4_976_067
DELETE_END = 4_976_790
BIN_BP = 1_000
DISPLAY_START = 4_803_501
DISPLAY_END = 5_144_387


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def quantile_match(source: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Map source ranks to the sorted empirical reference distribution."""
    source = np.asarray(source, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    ranks = rankdata(source, method="average") - 1.0
    q = ranks / max(source.size - 1, 1)
    ref_sorted = np.sort(reference)
    x = q * max(ref_sorted.size - 1, 0)
    lo = np.floor(x).astype(int)
    hi = np.ceil(x).astype(int)
    alpha = x - lo
    return (1.0 - alpha) * ref_sorted[lo] + alpha * ref_sorted[hi]


def write_fasta(path: Path, name: str, sequence: str, width: int = 80) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(f">{name}\n")
        for i in range(0, len(sequence), width):
            handle.write(sequence[i : i + width] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--control-dnase", required=True, type=Path)
    parser.add_argument("--experimental-npz", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    genome = Fasta(str(args.fasta), as_raw=True, sequence_always_upper=True)
    deletion_length = DELETE_END - DELETE_START
    wt = str(genome[CHROM][START:END]).upper()
    wt_extended = str(genome[CHROM][START : END + deletion_length]).upper()
    if len(wt) != END - START:
        raise RuntimeError(f"unexpected reference length: {len(wt)}")
    rel_start = DELETE_START - START
    rel_end = DELETE_END - START
    # Released predictors require a fixed-size input.  Construct the derivative
    # chromosome by joining the left flank to sequence after the right cut and
    # extending the right edge by the deletion length.  This is a true deletion,
    # not an N-mask or reference-coordinate substitution.
    deletion = wt_extended[:rel_start] + wt_extended[rel_end:]
    if len(deletion) != END - START:
        raise RuntimeError(f"unexpected derivative sequence length: {len(deletion)}")

    with pyBigWig.open(str(args.control_dnase)) as bw:
        control_bp = np.asarray(bw.values(CHROM, START, END, numpy=True), dtype=np.float32)
        control_bp_extended = np.asarray(
            bw.values(CHROM, START, END + deletion_length, numpy=True), dtype=np.float32
        )
    control_bp = np.nan_to_num(control_bp, nan=0.0, posinf=0.0, neginf=0.0)
    control_bp_extended = np.nan_to_num(
        control_bp_extended, nan=0.0, posinf=0.0, neginf=0.0
    )
    deletion_control_bp = np.concatenate(
        [control_bp_extended[:rel_start], control_bp_extended[rel_end:]]
    )
    if deletion_control_bp.size != END - START:
        raise RuntimeError("deletion-aligned control DNase has the wrong length")
    control_kb = control_bp.reshape(-1, BIN_BP).mean(axis=1)
    deletion_control_kb = deletion_control_bp.reshape(-1, BIN_BP).mean(axis=1)

    truth = np.load(args.experimental_npz)
    required = {"wt", "deletion", "wt_endpoint", "deletion_endpoint", "start", "end", "bin_bp"}
    if not required.issubset(truth.files):
        raise RuntimeError(f"experimental NPZ missing keys: {sorted(required - set(truth.files))}")
    if int(truth["start"]) != START or int(truth["end"]) != END or int(truth["bin_bp"]) != BIN_BP:
        raise RuntimeError("experimental truth geometry does not match the preregistered context")

    deletion_endpoint = np.asarray(truth["deletion_endpoint"], dtype=np.float64)
    surrogate_kb = quantile_match(np.log1p(deletion_endpoint), control_kb)
    surrogate_bp_reference = np.repeat(surrogate_kb.astype(np.float32), BIN_BP)
    # The target-assisted profile exists only on the registered 4-Mb reference
    # interval.  Its derivative-coordinate version shifts signal across the
    # deletion and uses the unchanged control-DNase tail solely to restore the
    # required fixed input length at the far-right boundary.
    surrogate_bp = np.concatenate(
        [
            surrogate_bp_reference[:rel_start],
            surrogate_bp_reference[rel_end:],
            control_bp_extended[-deletion_length:],
        ]
    ).astype(np.float32)
    surrogate_kb_derivative = surrogate_bp.reshape(-1, BIN_BP).mean(axis=1)

    derivative_to_reference = START + np.arange(END - START, dtype=np.int64)
    derivative_to_reference[rel_start:] += deletion_length

    wt_fasta = args.out_dir / "BHLHE40_WT_chr3_3Mb_7Mb.fa"
    deletion_fasta = args.out_dir / "BHLHE40_DEL723_chr3_3Mb_7Mb.fa"
    write_fasta(wt_fasta, "chr3_3000000_7000000_WT_hg38", wt)
    write_fasta(
        deletion_fasta,
        "chr3_3000000_7000000_DEL4976067_4976790_hg38",
        deletion,
    )

    crop0 = (DISPLAY_START - START) // BIN_BP
    crop1 = int(np.ceil((DISPLAY_END - START) / BIN_BP))
    out_npz = args.out_dir / "COMMON_INPUTS_AND_TRUTH.npz"
    np.savez_compressed(
        out_npz,
        wt_sequence=np.frombuffer(wt.encode("ascii"), dtype=np.uint8),
        deletion_sequence=np.frombuffer(deletion.encode("ascii"), dtype=np.uint8),
        control_dnase_bp=control_bp,
        control_dnase_1kb=control_kb.astype(np.float32),
        deletion_control_dnase_bp=deletion_control_bp.astype(np.float32),
        deletion_control_dnase_1kb=deletion_control_kb.astype(np.float32),
        target_assisted_surrogate_1kb=surrogate_kb.astype(np.float32),
        target_assisted_surrogate_bp_reference=surrogate_bp_reference,
        target_assisted_surrogate_bp_derivative=surrogate_bp,
        target_assisted_surrogate_1kb_derivative=surrogate_kb_derivative,
        derivative_to_reference_bp=derivative_to_reference,
        observed_wt_1kb=np.asarray(truth["wt"], dtype=np.float32),
        observed_deletion_1kb=np.asarray(truth["deletion"], dtype=np.float32),
        observed_wt_endpoint_1kb=np.asarray(truth["wt_endpoint"], dtype=np.float32),
        observed_deletion_endpoint_1kb=np.asarray(truth["deletion_endpoint"], dtype=np.float32),
        observed_wt_crop_1kb=np.asarray(truth["wt"][crop0:crop1, crop0:crop1], dtype=np.float32),
        observed_deletion_crop_1kb=np.asarray(truth["deletion"][crop0:crop1, crop0:crop1], dtype=np.float32),
        start=np.int64(START),
        end=np.int64(END),
        display_start=np.int64(DISPLAY_START),
        display_end=np.int64(DISPLAY_END),
        delete_start=np.int64(DELETE_START),
        delete_end=np.int64(DELETE_END),
        bin_bp=np.int64(BIN_BP),
    )

    metadata = {
        "status": "COMPLETE",
        "genome": "hg38",
        "chromosome": CHROM,
        "context_0based_half_open": [START, END],
        "display_0based_half_open": [DISPLAY_START, DISPLAY_END],
        "deletion_0based_half_open": [DELETE_START, DELETE_END],
        "deletion_length_bp": deletion_length,
        "deletion_status": (
            "inferred clean SpCas9 cut-to-cut deletion; exact clone junction unavailable"
        ),
        "deletion_window_construction": (
            "fixed-length derivative chromosome: delete the inferred interval, "
            "join the flanks, and append an equal-length downstream reference tail"
        ),
        "coordinate_mapping": (
            "derivative offsets before the junction map directly to hg38; offsets "
            "at/after the junction map to hg38 plus 723 bp"
        ),
        "strict_de_novo_accessibility": (
            "the same control K562 DNase observation is used for WT and deletion; "
            "for the deletion it is shifted with its attached DNA across the junction"
        ),
        "target_assisted_surrogate": {
            "status": "NOT_DE_NOVO",
            "source": "experimental deletion Hi-TrAC 1D endpoint profile",
            "transformation": "rank/quantile match log1p endpoint profile to empirical control-DNase 1-kb distribution",
            "control_dnase_sum_1kb": float(control_kb.sum()),
            "surrogate_sum_1kb_reference": float(surrogate_kb.sum()),
            "surrogate_sum_1kb_derivative": float(surrogate_kb_derivative.sum()),
            "control_dnase_mean_1kb": float(control_kb.mean()),
            "surrogate_mean_1kb": float(surrogate_kb.mean()),
            "spearman_surrogate_vs_deletion_hitrac_1d": float(
                spearmanr(surrogate_kb, deletion_endpoint).statistic
            ),
        },
        "inputs": {
            "reference_fasta": {"path": str(args.fasta), "sha256": sha256(args.fasta)},
            "control_dnase": {"path": str(args.control_dnase), "sha256": sha256(args.control_dnase)},
            "experimental_truth": {"path": str(args.experimental_npz), "sha256": sha256(args.experimental_npz)},
        },
        "outputs": {
            "common_npz": {"path": str(out_npz), "sha256": sha256(out_npz)},
            "wt_fasta": {"path": str(wt_fasta), "sha256": sha256(wt_fasta)},
            "deletion_fasta": {"path": str(deletion_fasta), "sha256": sha256(deletion_fasta)},
        },
    }
    (args.out_dir / "COMMON_INPUTS_PROVENANCE.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
