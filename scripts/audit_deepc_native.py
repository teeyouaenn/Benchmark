#!/usr/bin/env python3
"""Validate the official DeepC native center-pole text outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def locate(out_dir: Path, condition: str) -> Path:
    matches = sorted((out_dir / condition).glob("class_predictions_*.txt"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one {condition} output, found {matches}")
    return matches[0]


def read_native(path: Path) -> tuple[np.ndarray, np.ndarray]:
    coords = []
    values = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 204:
                raise RuntimeError(f"{path}: expected 204 fields, found {len(fields)}")
            coords.append([int(fields[1]), int(fields[2])])
            values.append([float(value) for value in fields[3:]])
    coordinates = np.asarray(coords, dtype=np.int64)
    predictions = np.asarray(values, dtype=np.float32)
    if predictions.ndim != 2 or predictions.shape[1] != 201:
        raise RuntimeError(f"unexpected DeepC shape: {predictions.shape}")
    if not np.isfinite(predictions).all():
        raise RuntimeError(f"non-finite values in {path}")
    return coordinates, predictions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    paths = {condition: locate(args.out_dir, condition) for condition in ("wt", "deletion")}
    parsed = {condition: read_native(path) for condition, path in paths.items()}
    npz_path = args.out_dir / "DEEPC_NATIVE_POLES.npz"
    np.savez_compressed(
        npz_path,
        wt_coordinates=parsed["wt"][0],
        wt_predictions=parsed["wt"][1],
        deletion_coordinates=parsed["deletion"][0],
        deletion_predictions=parsed["deletion"][1],
        bin_bp=np.int64(5_000),
        offsets_bp=np.arange(-500_000, 500_001, 5_000, dtype=np.int64),
    )
    audit = {
        "status": "COMPLETE",
        "model": "official DeepC K562 5-kb checkpoint",
        "input": "DNA only",
        "k562_matched": True,
        "native_scientific_scale": "normalized center-anchored contact profile",
        "native_bin_bp": 5_000,
        "native_pole_length": 201,
        "native_files": {
            condition: {
                "path": str(path),
                "sha256": sha256(path),
                "number_of_center_poles": int(parsed[condition][1].shape[0]),
                "minimum": float(parsed[condition][1].min()),
                "maximum": float(parsed[condition][1].max()),
                "mean": float(parsed[condition][1].mean()),
                "standard_deviation": float(parsed[condition][1].std()),
            }
            for condition, path in paths.items()
        },
        "compressed_native_copy_sha256": sha256(npz_path),
    }
    (args.out_dir / "DEEPC_RUN_AUDIT.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
