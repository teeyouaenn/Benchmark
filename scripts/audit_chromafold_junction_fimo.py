#!/usr/bin/env python3
"""Run exact FIMO on every CTCF-motif window crossing the new junction."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np


DELETE_START = 4_976_067
FLANK_BP = 100


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fimo", required=True, type=Path)
    parser.add_argument("--meme-motifs", required=True, type=Path)
    parser.add_argument("--common-inputs", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--out-audit", required=True, type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    with np.load(args.common_inputs, allow_pickle=False) as common:
        sequence = bytes(common["deletion_sequence"]).decode("ascii")
        common_start = int(common["start"])
    junction_offset = DELETE_START - common_start
    segment_start_offset = junction_offset - FLANK_BP
    segment_end_offset = junction_offset + FLANK_BP
    segment = sequence[segment_start_offset:segment_end_offset]
    if len(segment) != 2 * FLANK_BP:
        raise RuntimeError("junction segment has the wrong length")
    segment_start_bp = common_start + segment_start_offset
    fasta = args.out_dir / "BHLHE40_DERIVATIVE_JUNCTION_200BP.fa"
    fasta.write_text(
        f">chr3_derivative_{segment_start_bp}_{segment_start_bp + len(segment)}\n{segment}\n",
        encoding="ascii",
    )

    version = subprocess.check_output([str(args.fimo), "--version"], text=True).strip()
    command = [
        str(args.fimo),
        "--text",
        "--thresh",
        "1e-4",
        "--max-stored-scores",
        "1000000",
        str(args.meme_motifs),
        str(fasta),
    ]
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    table_path = args.out_dir / "FIMO_JUNCTION.tsv"
    table_path.write_text(completed.stdout, encoding="utf-8")

    data_lines = [line for line in completed.stdout.splitlines() if line and not line.startswith("#")]
    hits: list[dict[str, object]] = []
    if data_lines:
        reader = csv.DictReader(io.StringIO("\n".join(data_lines)), delimiter="\t")
        for row in reader:
            start_0based = segment_start_bp + int(row["start"]) - 1
            end_0based = segment_start_bp + int(row["stop"])
            hits.append(
                {
                    "motif_id": row["motif_id"],
                    "strand": row["strand"],
                    "p_value": float(row["p-value"]),
                    "score": float(row["score"]),
                    "derivative_interval_0based_half_open": [start_0based, end_0based],
                    "crosses_junction": start_0based < DELETE_START < end_0based,
                }
            )
    crossing = [hit for hit in hits if hit["crosses_junction"]]
    if crossing:
        raise RuntimeError(f"FIMO detected retained motif(s) crossing the junction: {crossing}")

    payload = {
        "status": "PASS",
        "verdict": "No FIMO hit at p <= 1e-4 crosses the inferred derivative junction.",
        "fimo_version": version,
        "threshold": 1e-4,
        "junction_derivative_0based": DELETE_START,
        "segment_0based_half_open": [segment_start_bp, segment_start_bp + len(segment)],
        "significant_hits_in_200bp_segment": hits,
        "significant_hits_crossing_junction": crossing,
        "hashes": {
            "meme_motifs_sha256": sha256(args.meme_motifs),
            "junction_fasta_sha256": sha256(fasta),
            "junction_sequence_sha256": sha256_bytes(segment.encode("ascii")),
            "fimo_table_sha256": sha256(table_path),
        },
    }
    args.out_audit.parent.mkdir(parents=True, exist_ok=True)
    args.out_audit.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
