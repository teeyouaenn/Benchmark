#!/usr/bin/env python3
"""Fail-closed audit of the fixed-length BHLHE40 derivative inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


CHROM = "chr3"
START = 3_000_000
END = 7_000_000
DELETE_START = 4_976_067
DELETE_END = 4_976_790
DELETE_LEN = DELETE_END - DELETE_START


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def max_abs(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(left, dtype=np.float64) - np.asarray(right, dtype=np.float64))))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--common-inputs", required=True, type=Path)
    parser.add_argument("--reference-fasta", type=Path)
    parser.add_argument("--control-dnase", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    with np.load(args.common_inputs, allow_pickle=False) as data:
        wt = bytes(data["wt_sequence"])
        deletion = bytes(data["deletion_sequence"])
        control = np.asarray(data["control_dnase_bp"], dtype=np.float32)
        deletion_control = np.asarray(data["deletion_control_dnase_bp"], dtype=np.float32)
        mapping = np.asarray(data["derivative_to_reference_bp"], dtype=np.int64)
        geometry = {
            "start": int(data["start"]),
            "end": int(data["end"]),
            "delete_start": int(data["delete_start"]),
            "delete_end": int(data["delete_end"]),
        }

    expected_geometry = {
        "start": START,
        "end": END,
        "delete_start": DELETE_START,
        "delete_end": DELETE_END,
    }
    if geometry != expected_geometry:
        raise RuntimeError(f"unexpected derivative geometry: {geometry}")
    if len(wt) != END - START or len(deletion) != END - START:
        raise RuntimeError("DNA arrays are not exactly 4 Mb")
    if control.shape != (END - START,) or deletion_control.shape != (END - START,):
        raise RuntimeError("DNase arrays are not exactly 4 Mb")

    left = DELETE_START - START
    right = DELETE_END - START
    overlap_end = (END - START) - DELETE_LEN
    checks = {
        "dna_prejunction_exact": deletion[:left] == wt[:left],
        "dna_downstream_overlap_exact": deletion[left:overlap_end] == wt[right:],
        "dnase_prejunction_max_abs": max_abs(deletion_control[:left], control[:left]),
        "dnase_downstream_overlap_max_abs": max_abs(
            deletion_control[left:overlap_end], control[right:]
        ),
        "mapping_prejunction_exact": bool(
            np.array_equal(mapping[:left], START + np.arange(left, dtype=np.int64))
        ),
        "mapping_downstream_exact": bool(
            np.array_equal(
                mapping[left:], START + np.arange(left, END - START, dtype=np.int64) + DELETE_LEN
            )
        ),
    }
    if not checks["dna_prejunction_exact"] or not checks["dna_downstream_overlap_exact"]:
        raise RuntimeError("derivative DNA does not implement the registered 723-bp deletion")
    if checks["dnase_prejunction_max_abs"] != 0.0 or checks["dnase_downstream_overlap_max_abs"] != 0.0:
        raise RuntimeError("derivative DNase is not shifted with its attached reference DNA")
    if not checks["mapping_prejunction_exact"] or not checks["mapping_downstream_exact"]:
        raise RuntimeError("derivative-to-reference mapping is inconsistent")

    source_checks: dict[str, object] = {}
    if args.reference_fasta is not None:
        from pyfaidx import Fasta

        genome = Fasta(str(args.reference_fasta), as_raw=True, sequence_always_upper=True)
        extended = str(genome[CHROM][START : END + DELETE_LEN]).upper()
        expected_wt = extended[: END - START]
        expected_deletion = extended[:left] + extended[right:]
        source_checks["reference_wt_exact"] = wt.decode("ascii") == expected_wt
        source_checks["reference_deletion_exact"] = deletion.decode("ascii") == expected_deletion
        if not all(source_checks.values()):
            raise RuntimeError("DNA arrays do not match the declared hg38 source")
    if args.control_dnase is not None:
        import pyBigWig

        with pyBigWig.open(str(args.control_dnase)) as bigwig:
            try:
                values = bigwig.values(CHROM, START, END + DELETE_LEN, numpy=True)
            except TypeError:
                values = bigwig.values(CHROM, START, END + DELETE_LEN)
            extended = np.asarray(values, dtype=np.float32)
        extended = np.nan_to_num(extended, nan=0.0, posinf=0.0, neginf=0.0)
        expected_control = extended[: END - START]
        expected_deletion_control = np.concatenate([extended[:left], extended[right:]])
        source_checks["control_dnase_max_abs"] = max_abs(control, expected_control)
        source_checks["deletion_dnase_max_abs"] = max_abs(
            deletion_control, expected_deletion_control
        )
        if source_checks["control_dnase_max_abs"] != 0.0 or source_checks["deletion_dnase_max_abs"] != 0.0:
            raise RuntimeError("DNase arrays do not match the declared K562 bigWig")

    payload = {
        "status": "PASS",
        "deletion_interpretation": (
            "clean 723-bp SpCas9 cut-to-cut derivative inferred from the two published guides; "
            "the exact clone-resolved repair junction is not published"
        ),
        "geometry_0based_half_open": geometry,
        "checks": checks,
        "source_checks": source_checks,
        "hashes": {
            "wt_sequence_sha256": digest_bytes(wt),
            "deletion_sequence_sha256": digest_bytes(deletion),
            "control_dnase_float32_sha256": digest_bytes(control.tobytes()),
            "deletion_control_dnase_float32_sha256": digest_bytes(deletion_control.tobytes()),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
