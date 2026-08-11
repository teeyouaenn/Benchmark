#!/usr/bin/env python3
"""Run official EPCOT HFF Micro-C head with K562 accessibility inputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch


COMMON_START = 3_000_000
INPUT_START = 4_674_000
INPUT_END = 5_274_000
OUTPUT_START = 4_724_000
OUTPUT_END = 5_224_000
BIN_BP = 1_000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def one_hot_bins(sequence: np.ndarray) -> np.ndarray:
    lookup = np.full(256, -1, dtype=np.int8)
    for index, base in enumerate(b"ACGT"):
        lookup[base] = index
    labels = lookup[sequence].reshape(600, 1_000)
    out = np.zeros((600, 4, 1_000), dtype=np.float32)
    for channel in range(4):
        out[:, channel, :] = labels == channel
    return out


def make_input(sequence: np.ndarray, accessibility: np.ndarray) -> np.ndarray:
    seq = one_hot_bins(sequence)
    acc = np.asarray(accessibility, dtype=np.float32).reshape(600, 1, 1_000)
    return np.concatenate([seq, acc], axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--pretrain-checkpoint", required=True, type=Path)
    parser.add_argument("--common-npz", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.repository / "COP"))

    from microc.model import build_pretrain_model_microc  # noqa: PLC0415

    config = SimpleNamespace(
        bins=600,
        crop=50,
        nheads=4,
        hidden_dim=512,
        embed_dim=256,
        dim_feedforward=1024,
        enc_layers=1,
        dec_layers=2,
        dropout=0.2,
        fine_tune=False,
        trunk="transformer",
        pretrain_path=str(args.pretrain_checkpoint),
    )
    device = torch.device(args.device)
    model = build_pretrain_model_microc(config)
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    missing, unexpected = model.load_state_dict(checkpoint, strict=False)
    if missing:
        raise RuntimeError(
            f"EPCOT checkpoint mismatch; missing={missing}, unexpected={unexpected}"
        )
    model.to(device).eval()

    common = np.load(args.common_npz)
    start = INPUT_START - COMMON_START
    end = INPUT_END - COMMON_START
    inputs = {
        "wt_control_dnase": make_input(
            np.asarray(common["wt_sequence"])[start:end],
            np.asarray(common["control_dnase_bp"])[start:end],
        ),
        "deletion_control_dnase": make_input(
            np.asarray(common["deletion_sequence"])[start:end],
            np.asarray(common["deletion_control_dnase_bp"])[start:end],
        ),
        "deletion_target_assisted_hitrac1d": make_input(
            np.asarray(common["deletion_sequence"])[start:end],
            np.asarray(common["target_assisted_surrogate_bp_derivative"])[start:end],
        ),
    }
    predictions = {}
    with torch.inference_mode():
        for condition, array in inputs.items():
            tensor = torch.from_numpy(array[None]).to(device)
            output = model(tensor).squeeze(-1).squeeze(0)
            prediction = output.detach().cpu().numpy().astype(np.float32)
            if prediction.shape != (125_250,):
                raise RuntimeError(f"unexpected EPCOT output: {prediction.shape}")
            if not np.isfinite(prediction).all():
                raise RuntimeError("non-finite EPCOT prediction")
            predictions[condition] = prediction

    native_npz = args.out_dir / "EPCOT_NATIVE_UPPER_TRIANGLE.npz"
    np.savez_compressed(
        native_npz,
        **predictions,
        start=np.int64(OUTPUT_START),
        end=np.int64(OUTPUT_END),
        bin_bp=np.int64(BIN_BP),
        matrix_bins=np.int64(500),
    )
    audit = {
        "status": "COMPLETE",
        "model": "official EPCOT HFF Micro-C 1-kb checkpoint",
        "input": "DNA plus accessibility",
        "k562_matched": "accessibility only; output head is HFF Micro-C",
        "strict_de_novo_arm": (
            "WT and deletion use the same measured K562 control DNase; signal is "
            "shifted with its attached DNA across the deletion junction"
        ),
        "target_assisted_arm": (
            "experimental deletion Hi-TrAC 1D quantile-matched to the control-DNase "
            "distribution; this arm is not de novo and is ineligible for model ranking"
        ),
        "native_scientific_scale": "observed/expected ratio",
        "input_interval_0based_half_open": [INPUT_START, INPUT_END],
        "output_interval_0based_half_open": [OUTPUT_START, OUTPUT_END],
        "native_bin_bp": BIN_BP,
        "native_geometry": "flattened upper triangle of a 500 by 500 symmetric map",
        "checkpoint_sha256": sha256(args.checkpoint),
        "pretrain_checkpoint_sha256": sha256(args.pretrain_checkpoint),
        "ignored_legacy_checkpoint_keys": sorted(unexpected),
        "common_input_sha256": sha256(args.common_npz),
        "native_npz_sha256": sha256(native_npz),
        "summary": {
            condition: {
                "minimum": float(value.min()),
                "maximum": float(value.max()),
                "mean": float(value.mean()),
                "standard_deviation": float(value.std()),
            }
            for condition, value in predictions.items()
        },
    }
    (args.out_dir / "EPCOT_RUN_AUDIT.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
