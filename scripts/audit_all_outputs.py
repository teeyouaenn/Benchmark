#!/usr/bin/env python3
"""Fail-closed audit of the seven native Figure 3D prediction products."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "results" / "native"
OUT = ROOT / "results"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def finite(name: str, array: np.ndarray) -> None:
    if not np.isfinite(array).all():
        raise RuntimeError(f"{name} contains non-finite values")


def audit_npz(path: Path, expected: dict[str, tuple[int, ...]]) -> dict:
    with np.load(path, allow_pickle=False) as data:
        for key, shape in expected.items():
            if key not in data:
                raise KeyError(f"{path}: missing {key}")
            if data[key].shape != shape:
                raise ValueError(f"{path}:{key}: {data[key].shape} != {shape}")
            finite(f"{path}:{key}", data[key])
        return {key: list(data[key].shape) for key in expected}


def main() -> None:
    rows: list[dict[str, object]] = []

    akita_path = NATIVE / "akita_v2" / "preds.h5"
    with h5py.File(akita_path, "r") as handle:
        if handle["preds"].shape != (2, 130305, 5):
            raise ValueError("Unexpected AkitaV2 prediction shape")
        finite("AkitaV2", handle["preds"][:])
    rows.append({
        "model": "AkitaV2", "primary_channel": "HFF", "input": "DNA",
        "native_bin_bp": 2048, "native_geometry": "512-bin square reconstructed from k>=2 upper triangle",
        "native_scale": "processed log observed/expected", "absolute_counts": False,
        "k562_matched": False, "primary_de_novo": True, "output_sha256": sha256(akita_path),
    })

    deepc_path = NATIVE / "deepc" / "DEEPC_NATIVE_POLES.npz"
    audit_npz(deepc_path, {"wt_predictions": (83, 201), "deletion_predictions": (83, 201)})
    rows.append({
        "model": "DeepC", "primary_channel": "K562", "input": "DNA",
        "native_bin_bp": 5000, "native_geometry": "83 center-anchored poles x 201 offsets",
        "native_scale": "normalized center-anchored interaction profile", "absolute_counts": False,
        "k562_matched": True, "primary_de_novo": True, "output_sha256": sha256(deepc_path),
    })

    orca_path = NATIVE / "orca" / "ORCA_NATIVE_MATRICES.npz"
    audit_npz(orca_path, {"hff_wt": (250, 250), "hff_deletion": (250, 250)})
    rows.append({
        "model": "Orca", "primary_channel": "HFF", "input": "DNA",
        "native_bin_bp": 4000, "native_geometry": "250 x 250 dense symmetric matrix",
        "native_scale": "normalized contact enrichment", "absolute_counts": False,
        "k562_matched": False, "primary_de_novo": True, "output_sha256": sha256(orca_path),
    })

    epcot_path = NATIVE / "epcot" / "EPCOT_NATIVE_UPPER_TRIANGLE.npz"
    audit_npz(epcot_path, {"wt_control_dnase": (125250,), "deletion_control_dnase": (125250,),
                           "deletion_target_assisted_hitrac1d": (125250,)})
    rows.append({
        "model": "EPCOT", "primary_channel": "HFF Micro-C head", "input": "DNA + K562 control DNase",
        "native_bin_bp": 1000, "native_geometry": "500 x 500 dense symmetric matrix",
        "native_scale": "model-native predicted O/E scale", "absolute_counts": False,
        "k562_matched": "accessibility only", "primary_de_novo": True, "output_sha256": sha256(epcot_path),
    })

    chromafold_path = NATIVE / "chromafold" / "CHROMAFOLD_NATIVE_VSTRIPES.npz"
    audit_npz(chromafold_path, {"wt_control_dnase": (36, 400), "deletion_control_dnase": (36, 400),
                                "deletion_hitrac1d_assisted": (36, 400)})
    rows.append({
        "model": "ChromaFold motif", "primary_channel": "motif/no-coaccess checkpoint",
        "input": "K562 DNase proxy + hg38 CTCF motif", "native_bin_bp": 10000,
        "native_geometry": "36 center V-stripes x 400 partner offsets", "native_scale": "HiC-DC+ Z-score",
        "absolute_counts": False, "k562_matched": False, "primary_de_novo": True,
        "output_sha256": sha256(chromafold_path),
    })

    alpha_path = NATIVE / "alphagenome" / "ALPHAGENOME_NATIVE_TRACKDATA.npz"
    audit_npz(alpha_path, {"reference_values": (512, 512, 28), "deletion_values": (512, 512, 28)})
    rows.append({
        "model": "AlphaGenome", "primary_channel": "HFFc6 Micro-C", "input": "DNA",
        "native_bin_bp": 2048, "native_geometry": "512 x 512 dense symmetric matrix",
        "native_scale": "relative contact frequency", "absolute_counts": False,
        "k562_matched": False, "primary_de_novo": True, "output_sha256": sha256(alpha_path),
    })

    chimaera_path = NATIVE / "chimaera" / "CHIMAERA_NATIVE_ROTATED_MAPS.npz"
    audit_npz(chimaera_path, {"wt": (4, 32, 128), "deletion": (4, 32, 128)})
    rows.append({
        "model": "Chimaera", "primary_channel": "human release", "input": "DNA",
        "native_bin_bp": 2048, "native_geometry": "4 overlapping 32 x 128 rotated distance images",
        "native_scale": "standardized log-distance residual", "absolute_counts": False,
        "k562_matched": False, "primary_de_novo": True, "output_sha256": sha256(chimaera_path),
    })

    if len(rows) != 7:
        raise AssertionError("Expected seven models")
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "MODEL_INPUT_OUTPUT_AUDIT.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "status": "PASS", "models_complete": 7, "rows": rows,
        "scientific_boundary": "Native outputs are normalized structural predictions, not absolute PET counts.",
        "target_assisted_boundary": "Hi-TrAC-1D-assisted EPCOT/ChromaFold outputs are sensitivity analyses and ineligible for de novo ranking.",
    }
    active_files = sorted(
        path
        for path in NATIVE.rglob("*")
        if path.is_file() and "superseded_center_sampled" not in path.parts
    )
    manifest_path = OUT / "NATIVE_OUTPUT_SHA256SUMS.txt"
    manifest_path.write_text(
        "".join(
            f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}\n"
            for path in active_files
        ),
        encoding="utf-8",
    )
    payload["active_native_artifacts_hashed"] = len(active_files)
    payload["sha256_manifest"] = str(manifest_path)
    (OUT / "MODEL_INPUT_OUTPUT_AUDIT.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({"status": "PASS", "models_complete": 7, "csv": str(csv_path)}, indent=2))


if __name__ == "__main__":
    main()
