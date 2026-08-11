#!/usr/bin/env python3
"""Run the released ChromaFold motif/no-coaccessibility checkpoint.

The released model expects pseudobulk scATAC counts plus a sequence-derived
CTCF motif-score track at 50-bp resolution.  For this cross-assay experiment,
K562 DNase is explicitly treated as an accessibility proxy.  The primary
deletion arm shifts the unchanged control-DNase signal with the derivative
chromosome.  Motif intervals, rather than already-binned motif values, are
transformed onto the derivative chromosome before exact 50-bp rasterization.
A target-assisted Hi-TrAC-1D surrogate is saved separately and is never labelled
de novo.
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


def sample_released_motif(motif_chr, reference_centers: np.ndarray) -> np.ndarray:
    indices = np.clip(
        reference_centers.astype(np.int64) // FINE_BP,
        0,
        motif_chr.shape[1] - 1,
    )
    selected = motif_chr[0, indices]
    if hasattr(selected, "toarray"):
        selected = selected.toarray()
    return np.asarray(selected).reshape(-1).astype(np.float32)


def transform_motif_hits(
    starts_1based: np.ndarray,
    ends_1based: np.ndarray,
    deletion: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map intact reference motifs onto the derivative chromosome.

    The deleted 0-based interval corresponds to 1-based bases
    DELETE_START+1 through DELETE_END.  Motifs overlapping any deleted base are
    disrupted and removed; intact downstream motifs shift left by DELETE_LEN.
    """
    starts = np.asarray(starts_1based, dtype=np.int64).copy()
    ends = np.asarray(ends_1based, dtype=np.int64).copy()
    keep = np.ones(starts.size, dtype=bool)
    if deletion:
        deleted_first_1based = DELETE_START + 1
        deleted_last_1based = DELETE_END
        overlaps_deleted = (ends >= deleted_first_1based) & (starts <= deleted_last_1based)
        keep &= ~overlaps_deleted
        downstream = starts > deleted_last_1based
        starts[downstream] -= DELETE_LEN
        ends[downstream] -= DELETE_LEN
    return starts, ends, keep


def rasterize_motif_hits(
    hits: dict[str, np.ndarray],
    first_center: int,
    last_center: int,
    deletion: bool,
) -> tuple[np.ndarray, dict[str, int]]:
    """Reproduce the released 50-bp motif feature on one coordinate system."""
    if first_center % FINE_BP != FINE_BP // 2 or last_center % FINE_BP != FINE_BP // 2:
        raise RuntimeError("motif centers are not aligned to the released 50-bp grid")
    first_bin = first_center // FINE_BP
    last_bin = last_center // FINE_BP
    size = last_bin - first_bin + 1
    forward = np.zeros(size, dtype=np.float64)
    reverse = np.zeros(size, dtype=np.float64)
    starts, ends, keep = transform_motif_hits(
        hits["start_1based"], hits["end_1based"], deletion
    )
    strands = np.asarray(hits["strand_code"], dtype=np.int8)
    scores = np.asarray(hits["score"], dtype=np.float64)
    for strand_code, target in ((1, forward), (-1, reverse)):
        selected = np.flatnonzero(keep & (strands == strand_code))
        for index in selected:
            start = int(starts[index])
            end = int(ends[index])
            score = float(scores[index])
            candidate_first = max(first_bin, (start - FINE_BP) // FINE_BP)
            candidate_last = min(last_bin, end // FINE_BP)
            for bin_index in range(candidate_first, candidate_last + 1):
                overlap = max(
                    0,
                    min(end, bin_index * FINE_BP + FINE_BP)
                    - max(start, bin_index * FINE_BP)
                    + 1,
                )
                if overlap >= 10:
                    target[bin_index - first_bin] = score
    feature = np.log1p(np.where(reverse > 0, reverse, forward)).astype(np.float32)
    return feature, {
        "motifs_input": int(starts.size),
        "motifs_disrupted_by_deletion": int(np.count_nonzero(~keep)),
        "motifs_shifted_downstream": int(
            np.count_nonzero(
                keep & (np.asarray(hits["start_1based"], dtype=np.int64) > DELETE_END)
            )
        ) if deletion else 0,
        "nonzero_50bp_bins": int(np.count_nonzero(feature)),
    }


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
    motif_track: np.ndarray,
    motif_first_center: int,
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
        motif_offset = int(derivative_centers[0] - motif_first_center) // FINE_BP
        motif = motif_track[motif_offset : motif_offset + n_fine]
        if motif.size != n_fine:
            raise RuntimeError("motif track does not cover a requested ChromaFold window")
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
    parser.add_argument("--motif-hits", required=True, type=Path)
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
    with np.load(args.motif_hits, allow_pickle=False) as motif_resource:
        motif_hits = {
            key: np.asarray(motif_resource[key]).copy()
            for key in ("start_1based", "end_1based", "strand_code", "score")
        }
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
        n_fine = INPUT_BP // FINE_BP
        motif_first_center = int(centers.min()) - 2_000_000 + FINE_BP // 2
        motif_last_center = (
            int(centers.max()) - 2_000_000 + FINE_BP // 2 + (n_fine - 1) * FINE_BP
        )
        motif_reference, reference_motif_summary = rasterize_motif_hits(
            motif_hits, motif_first_center, motif_last_center, deletion=False
        )
        motif_deletion, deletion_motif_summary = rasterize_motif_hits(
            motif_hits, motif_first_center, motif_last_center, deletion=True
        )
        motif_centers = np.arange(
            motif_first_center, motif_last_center + 1, FINE_BP, dtype=np.int64
        )
        released_reference = sample_released_motif(motif_chr, motif_centers)
        reference_difference = np.abs(motif_reference - released_reference)
        if float(reference_difference.max()) > 1e-6:
            raise RuntimeError(
                "motif-level reconstruction does not reproduce the released WT feature: "
                f"max_abs={float(reference_difference.max())}"
            )
        legacy_deletion = sample_released_motif(
            motif_chr, derivative_to_reference(motif_centers)
        )
        legacy_difference = np.abs(motif_deletion - legacy_deletion)
        signal_start = int(centers.min()) - 2_000_000
        signal_end = int(centers.max()) + 2_010_000 + DELETE_LEN
        try:
            values = bw.values(CHROM, signal_start, signal_end, numpy=True)
        except TypeError:
            values = bw.values(CHROM, signal_start, signal_end)
        signal_bp = np.asarray(values, dtype=np.float64)
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
            condition_motif = (
                motif_reference if condition == "wt_control_dnase" else motif_deletion
            )
            x = make_inputs(
                centers,
                condition,
                signal_prefix,
                signal_start,
                condition_motif,
                motif_first_center,
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
                "intact motif intervals and unchanged control DNase transformed onto the "
                "inferred 723-bp cut-to-cut derivative chromosome"
            ),
            "deletion_hitrac1d_assisted": (
                "NOT DE NOVO: experimental deletion Hi-TrAC 1D rank/quantile "
                "matched to control DNase units"
            ),
        },
        "adaptations": [
            "bulk K562 DNase substitutes for the checkpoint's pseudobulk scATAC input",
            "the official whole-genome library-size normalization and log1p transform are retained",
            (
                "official AH104727 motif intervals are transformed before 50-bp rasterization; "
                "motifs disrupted by the deletion are removed"
            ),
            (
                "all motif windows crossing the new junction were audited against the three "
                "JASPAR 2022 CTCF PWMs and none approached the official retention boundary"
            ),
            "no scATAC coaccessibility is used, matching the selected no-coaccessibility checkpoint",
        ],
        "deletion_status": (
            "inferred clean SpCas9 cut-to-cut deletion; exact clone junction unavailable"
        ),
        "motif_derivative_audit": {
            "reference_reconstruction_bins": int(motif_reference.size),
            "reference_reconstruction_max_abs_difference": float(reference_difference.max()),
            "reference_reconstruction_mismatched_bins_at_1e-6": int(
                np.count_nonzero(reference_difference > 1e-6)
            ),
            "legacy_center_sampled_deletion_bins_changed": int(
                np.count_nonzero(legacy_difference > 1e-7)
            ),
            "legacy_center_sampled_deletion_max_abs_difference": float(
                legacy_difference.max()
            ),
            "reference": reference_motif_summary,
            "deletion": deletion_motif_summary,
        },
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
            "motif_hits_sha256": sha256(args.motif_hits),
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
