#!/usr/bin/env python3
"""Compare the corrected ChromaFold derivative against its superseded run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


CONDITIONS = (
    "wt_control_dnase",
    "deletion_control_dnase",
    "deletion_hitrac1d_assisted",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compare(left: np.ndarray, right: np.ndarray) -> dict[str, object]:
    difference = np.abs(left - right)
    return {
        "pearson": float(np.corrcoef(left.ravel(), right.ravel())[0, 1]),
        "mean_absolute_difference": float(difference.mean()),
        "max_absolute_difference": float(difference.max()),
        "exactly_equal": bool(np.array_equal(left, right)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--superseded", required=True, type=Path)
    parser.add_argument("--corrected", required=True, type=Path)
    parser.add_argument("--repeat", required=True, type=Path)
    parser.add_argument("--superseded-audit", required=True, type=Path)
    parser.add_argument("--corrected-audit", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    old = np.load(args.superseded, allow_pickle=False)
    new = np.load(args.corrected, allow_pickle=False)
    repeat = np.load(args.repeat, allow_pickle=False)
    old_audit = json.loads(args.superseded_audit.read_text(encoding="utf-8"))
    new_audit = json.loads(args.corrected_audit.read_text(encoding="utf-8"))

    reproducibility = {
        condition: compare(new[condition], repeat[condition]) for condition in CONDITIONS
    }
    if not all(item["exactly_equal"] for item in reproducibility.values()):
        raise RuntimeError("corrected ChromaFold inference is not exactly reproducible")
    if old_audit["input_summaries"]["wt_control_dnase"] != new_audit["input_summaries"][
        "wt_control_dnase"
    ]:
        raise RuntimeError("WT input changed during the deletion-specific correction")

    payload = {
        "status": "PASS",
        "verdict": (
            "The prior center-sampled 50-bp deletion motif track is superseded. The corrected "
            "run transforms official AH104727 motif intervals onto the derivative chromosome "
            "before rasterization; WT inputs are unchanged and corrected inference is exactly reproducible."
        ),
        "superseded_sha256": sha256(args.superseded),
        "corrected_sha256": sha256(args.corrected),
        "comparison": {
            condition: compare(old[condition], new[condition]) for condition in CONDITIONS
        },
        "repeat_reproducibility": reproducibility,
        "input_invariants": {
            "wt_input_summary_exact": True,
            "deletion_accessibility_mean_exact": (
                old_audit["input_summaries"]["deletion_control_dnase"]["accessibility_mean"]
                == new_audit["input_summaries"]["deletion_control_dnase"]["accessibility_mean"]
            ),
            "deletion_accessibility_sd_exact": (
                old_audit["input_summaries"]["deletion_control_dnase"]["accessibility_sd"]
                == new_audit["input_summaries"]["deletion_control_dnase"]["accessibility_sd"]
            ),
        },
        "motif_derivative_audit": new_audit["motif_derivative_audit"],
    }
    if not all(payload["input_invariants"].values()):
        raise RuntimeError("a non-motif input changed during correction")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
