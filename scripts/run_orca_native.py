#!/usr/bin/env python3
"""Run official Orca 1-Mb H1-ESC and HFF checkpoints on WT/deletion DNA."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch


COMMON_START = 3_000_000
WINDOW_START = 4_474_000
WINDOW_END = 5_474_000
BIN_BP = 4_000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def one_hot(sequence: np.ndarray) -> np.ndarray:
    lookup = np.full(256, -1, dtype=np.int8)
    for index, base in enumerate(b"ACGT"):
        lookup[base] = index
    labels = lookup[sequence]
    out = np.zeros((sequence.size, 4), dtype=np.float32)
    valid = labels >= 0
    out[np.nonzero(valid)[0], labels[valid]] = 1.0
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--common-npz", required=True, type=Path)
    parser.add_argument("--orca-path", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--device", default="cuda:0")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.orca_path))
    from orca_models import H1esc_1M, Hff_1M  # noqa: PLC0415

    common = np.load(args.common_npz)
    start = WINDOW_START - COMMON_START
    end = WINDOW_END - COMMON_START
    sequences = {
        "wt": np.asarray(common["wt_sequence"])[start:end],
        "deletion": np.asarray(common["deletion_sequence"])[start:end],
    }
    if any(sequence.size != 1_000_000 for sequence in sequences.values()):
        raise RuntimeError("Orca input must contain exactly 1,000,000 bases")

    device = torch.device(args.device)
    models = {"h1esc": H1esc_1M(), "hff": Hff_1M()}
    outputs: dict[str, dict[str, np.ndarray]] = {}
    with torch.inference_mode():
        for model_name, model in models.items():
            model = model.to(device).eval()
            outputs[model_name] = {}
            for condition, sequence in sequences.items():
                tensor = torch.from_numpy(one_hot(sequence).T[None]).to(device)
                prediction = model(tensor).squeeze().detach().cpu().numpy().astype(np.float32)
                if prediction.shape != (250, 250):
                    raise RuntimeError(f"unexpected Orca output shape: {prediction.shape}")
                if not np.isfinite(prediction).all():
                    raise RuntimeError("non-finite Orca prediction")
                outputs[model_name][condition] = prediction
            model.to("cpu")
            del model
            if device.type == "cuda":
                torch.cuda.empty_cache()

    native_path = args.out_dir / "ORCA_NATIVE_OUTPUT.pth"
    torch.save(
        {
            "predictions": outputs,
            "chromosome": "chr3",
            "start": WINDOW_START,
            "end": WINDOW_END,
            "bin_bp": BIN_BP,
            "conditions": {
                "wt": "hg38 reference",
                "deletion": "inferred 723-bp clean deletion",
            },
        },
        native_path,
    )
    npz_path = args.out_dir / "ORCA_NATIVE_MATRICES.npz"
    np.savez_compressed(
        npz_path,
        h1esc_wt=outputs["h1esc"]["wt"],
        h1esc_deletion=outputs["h1esc"]["deletion"],
        hff_wt=outputs["hff"]["wt"],
        hff_deletion=outputs["hff"]["deletion"],
        start=np.int64(WINDOW_START),
        end=np.int64(WINDOW_END),
        bin_bp=np.int64(BIN_BP),
    )
    audit = {
        "status": "COMPLETE",
        "model": "official Orca 1-Mb checkpoints",
        "cell_types": ["H1-ESC", "HFF"],
        "k562_matched": False,
        "input": "DNA only",
        "native_scientific_scale": "normalized contact enrichment",
        "input_window_0based_half_open": [WINDOW_START, WINDOW_END],
        "native_bin_bp": BIN_BP,
        "native_shape": [250, 250],
        "common_input_sha256": sha256(args.common_npz),
        "native_pth_sha256": sha256(native_path),
        "native_npz_sha256": sha256(npz_path),
        "summary": {
            name: {
                condition: {
                    "minimum": float(array.min()),
                    "maximum": float(array.max()),
                    "mean": float(array.mean()),
                    "standard_deviation": float(array.std()),
                }
                for condition, array in condition_outputs.items()
            }
            for name, condition_outputs in outputs.items()
        },
    }
    (args.out_dir / "ORCA_RUN_AUDIT.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
