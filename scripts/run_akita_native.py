#!/usr/bin/env python3
"""Run official AkitaV2 fold-2 human checkpoint on WT/deletion DNA."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("TF_FORCE_GPU_ALLOW_GROWTH", "true")
os.environ.setdefault("TF_USE_LEGACY_KERAS", "1")

import h5py
import numpy as np


COMMON_START = 3_000_000
INPUT_START = 4_318_640
INPUT_END = 5_629_360
BIN_BP = 2_048
OUTPUT_START = INPUT_START + 64 * BIN_BP
OUTPUT_END = INPUT_END - 64 * BIN_BP


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
    parser.add_argument("--repository", required=True, type=Path)
    parser.add_argument("--params", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--common-npz", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.repository))

    from basenji import seqnn  # noqa: PLC0415

    params = json.loads(args.params.read_text(encoding="utf-8"))
    model = seqnn.SeqNN(params["model"])
    model.restore(str(args.checkpoint))
    # Reverse-complement averaging is an official Akita inference option and
    # produces one prediction on the native scientific scale.
    model.build_ensemble(ensemble_rc=True, ensemble_shifts=[0])

    common = np.load(args.common_npz)
    start = INPUT_START - COMMON_START
    end = INPUT_END - COMMON_START
    input_batch = np.stack(
        [
            one_hot(np.asarray(common["wt_sequence"])[start:end]),
            one_hot(np.asarray(common["deletion_sequence"])[start:end]),
        ],
        axis=0,
    )
    if input_batch.shape != (2, 1_310_720, 4):
        raise RuntimeError(f"unexpected Akita input shape: {input_batch.shape}")
    prediction = np.asarray(model(input_batch), dtype=np.float32)
    if prediction.ndim != 3 or prediction.shape[0] != 2 or prediction.shape[2] != 5:
        raise RuntimeError(f"unexpected Akita output shape: {prediction.shape}")
    if not np.isfinite(prediction).all():
        raise RuntimeError("non-finite Akita prediction")

    native_h5 = args.out_dir / "preds.h5"
    with h5py.File(native_h5, "w") as handle:
        dataset = handle.create_dataset("preds", data=prediction, compression="gzip")
        dataset.attrs["condition_order"] = np.asarray(["wt", "deletion"], dtype="S")
        dataset.attrs["native_bin_bp"] = BIN_BP
        dataset.attrs["output_start_0based"] = OUTPUT_START
        dataset.attrs["output_end_0based"] = OUTPUT_END
    target_copy = args.out_dir / "targets.txt"
    target_copy.write_bytes(args.targets.read_bytes())
    audit = {
        "status": "COMPLETE",
        "model": "official AkitaV2 fold-2 human checkpoint",
        "input": "DNA only",
        "k562_matched": False,
        "target_channels": ["HFF", "H1-ESC", "GM12878", "IMR-90", "HCT116"],
        "inference_ensemble": "forward plus reverse-complement average, shift 0",
        "native_scientific_scale": "processed log observed/expected",
        "input_interval_0based_half_open": [INPUT_START, INPUT_END],
        "output_interval_0based_half_open": [OUTPUT_START, OUTPUT_END],
        "native_bin_bp": BIN_BP,
        "native_shape": list(prediction.shape),
        "params_sha256": sha256(args.params),
        "checkpoint_sha256": sha256(args.checkpoint),
        "common_input_sha256": sha256(args.common_npz),
        "native_hdf5_sha256": sha256(native_h5),
        "summary": {
            condition: {
                "minimum": float(prediction[index].min()),
                "maximum": float(prediction[index].max()),
                "mean": float(prediction[index].mean()),
                "standard_deviation": float(prediction[index].std()),
            }
            for index, condition in enumerate(("wt", "deletion"))
        },
    }
    (args.out_dir / "AKITAV2_RUN_AUDIT.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
