#!/usr/bin/env python3
"""Run the released ChromaFold motif/no-coaccessibility checkpoint.

The released model expects pseudobulk scATAC counts plus a sequence-derived
CTCF motif-score track at 50-bp resolution.  For this cross-assay experiment,
K562 DNase is explicitly treated as an accessibility proxy.  The primary
deletion arm shifts the unchanged control-DNase signal with the derivative
chromosome; a target-assisted Hi-TrAC-1D surrogate is saved separately and is
never labelled de novo.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyBigWig
import torch


CHROM = "chr3"
DISPLAY_START = 4_803_501
DISPLAY_END = 5_144_387
DELETE_START = 4_976_067
DELETE_END = 4_976_790
DELETE_LEN = DELETE_END - DELETE_START
COMMON_START = 3_000_000
COMMON_END = 7_000_000
INPUT_BP = 4_010_000
FINE_BP = 50
OUTPUT_BP = 10_000
EFFECTIVE_GENOME_BP = 2_913_022_398


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_head(path: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()


def derivative_to_reference(position: np.ndarray) -> np.ndarray:
    position = np.asarray(position, dtype=np.int64)
    return np.where(position < DELETE_START, position, position + DELETE_LEN)


def bw_means_from_buffer(
    reference_centers: np.ndarray,
    signal_prefix: np.ndarray,
    signal_start: int,
) -> np.ndarray:
    """Return exact 50-bp means from one preloaded base-resolution interval."""
    starts = reference_centers.astype(np.int64) - FINE_BP // 2 - signal_start
    ends = starts + FINE_BP
    if starts.min() < 0 or ends.max() >= signal_prefix.size:
        raise RuntimeError("reference-coordinate request exceeds preloaded DNase buffer")
    return ((signal_prefix[ends] - signal_prefix[starts]) / FINE_BP).astype(np.float32)


def sample_motif(motif_chr, reference_centers: np.ndarray) -> np.ndarray:
    indices = np.clip(
        reference_centers.astype(np.int64) // FINE_BP,
        0,
        motif_chr.shape[1] - 1,
    )
    selected = motif_chr[0, indices]
    if hasattr(selected, "toarray"):
        selected = selected.toarray()
    return np.asarray(selected).reshape(-1).astype(np.float32)


def surrogate_signal(
    derivative_centers: np.ndarray,
    surrogate_bp: np.ndarray,
    fallback: np.ndarray,
) -> np.ndarray:
    out = fallback.copy()
    inside = (derivative_centers >= COMMON_START) & (derivative_centers < COMMON_END)
    indices = derivative_centers[inside] - COMMON_START
    out[inside] = surrogate_bp[indices]
    return out


def make_inputs(
    centers: np.ndarray,
    condition: str,
    signal_prefix: np.ndarray,
    signal_start: int,
    motif_chr,
    dnase_scale: float,
    surrogate_bp: np.ndarray,
) -> torch.Tensor:
    inputs: list[np.ndarray] = []
    n_fine = INPUT_BP // FINE_BP
    for center in centers:
        start = int(center) - 2_000_000
        derivative_centers = start + FINE_BP // 2 + np.arange(n_fine) * FINE_BP
        if condition == "wt_control_dnase":
            reference_centers = derivative_centers
        else:
            reference_centers = derivative_to_reference(derivative_centers)

        dnase = bw_means_from_buffer(reference_centers, signal_prefix, signal_start)
        if condition == "deletion_hitrac1d_assisted":
            dnase = surrogate_signal(derivative_centers, surrogate_bp, dnase)
        accessibility = np.log1p(dnase / dnase_scale).astype(np.float32)
        motif = sample_motif(motif_chr, reference_centers)
        inputs.append(np.stack([accessibility, motif], axis=0))

    tensor = torch.from_numpy(np.stack(inputs, axis=0))
    if tensor.shape[1:] != (2, n_fine):
        raise RuntimeError(f"unexpected ChromaFold input shape {tuple(tensor.shape)}")
    return tensor


def infer(model: torch.nn.Module, inputs: torch.Tensor, device: torch.device) -> np.ndarray:
    outputs: list[np.ndarray] = []
    with torch.inference_mode():
        for begin in range(0, inputs.shape[0], 8):
            x = inputs[begin : begin + 8].to(device)
            left = model(x)
            right = torch.flip(model(torch.flip(x, dims=[2])), dims=[1])
            stripe = torch.cat([left, right], dim=1)
            outputs.append(stripe.cpu().numpy().astype(np.float32))
    result = np.concatenate(outputs, axis=0)
    if result.shape != (inputs.shape[0], 400):
        raise RuntimeError(f"unexpected ChromaFold output shape {result.shape}")
    if not np.isfinite(result).all():
        raise RuntimeError("non-finite ChromaFold predictions")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--motif-pickle", required=True, type=Path)
    parser.add_argument("--control-dnase", required=True, type=Path)
    parser.add_argument("--common-inputs", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--device", default="cuda:1")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(args.source / "chromafold"))
    from model_bulk_only import branch_pbulk  # noqa: PLC0415

    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    model = branch_pbulk().to(device)
    state = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    state = {key.removeprefix("module."): value for key, value in state.items()}
    incompatible = model.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise RuntimeError(f"checkpoint mismatch: {incompatible}")
    model.eval()

    with args.motif_pickle.open("rb") as handle:
        motif = pickle.load(handle)
    motif_chr = motif[CHROM]
    common = np.load(args.common_inputs)
    surrogate_bp = np.asarray(
        common["target_assisted_surrogate_bp_derivative"], dtype=np.float32
    )

    # Use the released ChromaFold whole-genome library-size convention.  For a
    # continuous bigWig, sumData/50 is the equivalent sum of 50-bp bin means.
    with pyBigWig.open(str(args.control_dnase)) as bw:
        header = bw.header()
        equivalent_libsize = float(header["sumData"]) / FINE_BP
        dnase_scale = equivalent_libsize * 150.0 / EFFECTIVE_GENOME_BP
        if not np.isfinite(dnase_scale) or dnase_scale <= 0:
            raise RuntimeError(f"invalid DNase scale {dnase_scale}")

        first_center = (DISPLAY_START // OUTPUT_BP) * OUTPUT_BP
        last_center = int(np.ceil(DISPLAY_END / OUTPUT_BP)) * OUTPUT_BP
        centers = np.arange(first_center, last_center + 1, OUTPUT_BP, dtype=np.int64)
        signal_start = int(centers.min()) - 2_000_000
        signal_end = int(centers.max()) + 2_010_000 + DELETE_LEN
        signal_bp = np.asarray(
            bw.values(CHROM, signal_start, signal_end, numpy=True), dtype=np.float64
        )
        signal_bp = np.nan_to_num(signal_bp, nan=0.0, posinf=0.0, neginf=0.0)
        signal_bp = np.maximum(signal_bp, 0.0)
        signal_prefix = np.concatenate([[0.0], np.cumsum(signal_bp, dtype=np.float64)])
        predictions = {}
        input_summaries = {}
        for condition in (
            "wt_control_dnase",
            "deletion_control_dnase",
            "deletion_hitrac1d_assisted",
        ):
            x = make_inputs(
                centers,
                condition,
                signal_prefix,
                signal_start,
                motif_chr,
                dnase_scale,
                surrogate_bp,
            )
            input_summaries[condition] = {
                "shape": list(x.shape),
                "accessibility_mean": float(x[:, 0].mean()),
                "accessibility_sd": float(x[:, 0].std()),
                "motif_mean": float(x[:, 1].mean()),
                "motif_sd": float(x[:, 1].std()),
            }
            predictions[condition] = infer(model, x, device)

    partner_offsets = np.concatenate(
        [np.arange(-200, 0), np.arange(1, 201)]
    ).astype(np.int64) * OUTPUT_BP
    output = args.out_dir / "CHROMAFOLD_NATIVE_VSTRIPES.npz"
    np.savez_compressed(
        output,
        center_positions_derivative=centers,
        center_positions_reference_wt=centers,
        center_positions_reference_deletion=derivative_to_reference(centers),
        partner_offsets_bp=partner_offsets,
        **predictions,
    )

    summary = {
        "status": "COMPLETE",
        "scientific_status": "ADAPTED_CROSS_ASSAY_NOT_NATIVE_K562_SCATAC",
        "checkpoint": {
            "name": "chromafold_CTCFmotif_noCoaccessibility.pth.tar",
            "sha256": sha256(args.checkpoint),
        },
        "source": {"git_head": git_head(args.source)},
        "native_prediction": {
            "description": "HiC-DC+ normalized Z-score V-stripe",
            "bin_bp": OUTPUT_BP,
            "input_span_bp": INPUT_BP,
            "output_geometry": "center by 400 partners; 200 left and 200 right",
            "shape_per_condition": list(next(iter(predictions.values())).shape),
            "output_file": str(output),
            "output_sha256": sha256(output),
        },
        "conditions": {
            "wt_control_dnase": "reference motif plus K562 control DNase proxy",
            "deletion_control_dnase": (
                "motif and the unchanged control DNase shifted with the inferred "
                "723-bp derivative chromosome"
            ),
            "deletion_hitrac1d_assisted": (
                "NOT DE NOVO: experimental deletion Hi-TrAC 1D rank/quantile "
                "matched to control DNase units"
            ),
        },
        "adaptations": [
            "bulk K562 DNase substitutes for the checkpoint's pseudobulk scATAC input",
            "the official whole-genome library-size normalization and log1p transform are retained",
            "reference CTCF motif scores are shifted with the deletion; the new junction is not rescanned",
            "no scATAC coaccessibility is used, matching the selected no-coaccessibility checkpoint",
        ],
        "deletion_status": (
            "inferred clean sgRNA-bounded deletion; exact clone junction unavailable"
        ),
        "dnase_equivalent_50bp_libsize": equivalent_libsize,
        "dnase_scale_factor": dnase_scale,
        "input_summaries": input_summaries,
        "prediction_summaries": {
            condition: {
                "min": float(array.min()),
                "max": float(array.max()),
                "mean": float(array.mean()),
                "sd": float(array.std()),
            }
            for condition, array in predictions.items()
        },
        "inputs": {
            "motif_pickle_sha256": sha256(args.motif_pickle),
            "control_dnase_sha256": sha256(args.control_dnase),
            "common_inputs_sha256": sha256(args.common_inputs),
        },
    }
    (args.out_dir / "CHROMAFOLD_RUN_AUDIT.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary["prediction_summaries"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
