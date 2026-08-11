#!/usr/bin/env python3
"""Render every released Figure 3D prediction as a square contact matrix.

Only geometry is standardized.  Predictor values remain in their released
scientific scales.  DeepC center poles and ChromaFold V-stripes are placed by
genomic coordinate and duplicate estimates are averaged.  Chimaera's native
distance-position image is mapped back to endpoint coordinates and overlapping
fragments are averaged.  Missing or unsupported cells remain NaN.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import fontManager


ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "results" / "native"
INPUTS = ROOT / "results" / "inputs" / "COMMON_INPUTS_AND_TRUTH.npz"
OUT_DIR = ROOT / "results" / "standardized"
FIGURES = ROOT / "figures"

DISPLAY_START = 4_803_501
DISPLAY_END = 5_144_387
EXPERIMENTAL_GRID_START = 4_803_000
EXPERIMENTAL_GRID_END = 5_145_000

PAPER_VMIN = 0.1
PAPER_VMAX = 0.8
PAPER_REDS = [
    "#ffffff", "#f4d2d2", "#f6b4b4", "#f79696", "#f97878",
    "#fa5959", "#fc3b3b", "#fe1d1d", "#ff0000",
]


@dataclass
class SquarePanel:
    slug: str
    name: str
    subtitle: str
    scale: str
    wt: np.ndarray
    deletion: np.ndarray
    start_bp: int
    end_bp: int
    bin_bp: int
    reconstruction: str
    source_geometry: str
    fixed_limits: tuple[float, float] | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_arial() -> None:
    arial = Path(r"C:\Windows\Fonts\arial.ttf")
    if arial.exists():
        fontManager.addfont(str(arial))
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )


def square_from_upper(values: np.ndarray, n: int, diagonal: int) -> np.ndarray:
    matrix = np.full((n, n), np.nan, dtype=np.float64)
    indices = np.triu_indices(n, k=diagonal)
    if len(indices[0]) != len(values):
        raise ValueError(f"Cannot place {len(values)} values into n={n}, k={diagonal}")
    matrix[indices] = values
    matrix[(indices[1], indices[0])] = values
    return matrix


def crop_square(
    matrix: np.ndarray, start_bp: int, bin_bp: int
) -> tuple[np.ndarray, int, int]:
    centers = start_bp + (np.arange(matrix.shape[0]) + 0.5) * bin_bp
    keep = (centers >= DISPLAY_START) & (centers < DISPLAY_END)
    indices = np.flatnonzero(keep)
    if len(indices) < 2:
        raise ValueError("Requested locus is outside model output")
    cropped = matrix[np.ix_(indices, indices)]
    crop_start = start_bp + int(indices[0]) * bin_bp
    crop_end = start_bp + (int(indices[-1]) + 1) * bin_bp
    return cropped, crop_start, crop_end


def lattice_positions(
    coordinates: np.ndarray, bin_bp: int
) -> np.ndarray:
    values = np.asarray(coordinates, dtype=np.int64)
    residues, counts = np.unique(np.mod(values, bin_bp), return_counts=True)
    residue = int(residues[np.argmax(counts)])
    first = int(np.ceil((DISPLAY_START - residue) / bin_bp) * bin_bp + residue)
    last = int(np.floor((DISPLAY_END - 1 - residue) / bin_bp) * bin_bp + residue)
    if last < first:
        raise ValueError("No model lattice positions overlap the display locus")
    return np.arange(first, last + bin_bp, bin_bp, dtype=np.int64)


def aggregate_pairs(
    left: np.ndarray,
    right: np.ndarray,
    values: np.ndarray,
    positions: np.ndarray,
) -> tuple[np.ndarray, dict[str, float | int]]:
    sums: dict[tuple[int, int], float] = defaultdict(float)
    counts: dict[tuple[int, int], int] = defaultdict(int)
    position_to_index = {int(value): index for index, value in enumerate(positions)}
    used = 0
    for a, b, value in zip(left, right, values, strict=True):
        a_int = int(a)
        b_int = int(b)
        if a_int not in position_to_index or b_int not in position_to_index:
            continue
        key = (min(a_int, b_int), max(a_int, b_int))
        sums[key] += float(value)
        counts[key] += 1
        used += 1

    matrix = np.full((len(positions), len(positions)), np.nan, dtype=np.float64)
    duplicate_spreads: list[float] = []
    grouped_values: dict[tuple[int, int], list[float]] = defaultdict(list)
    for a, b, value in zip(left, right, values, strict=True):
        a_int = int(a)
        b_int = int(b)
        if a_int in position_to_index and b_int in position_to_index:
            grouped_values[(min(a_int, b_int), max(a_int, b_int))].append(float(value))

    for key, total in sums.items():
        i = position_to_index[key[0]]
        j = position_to_index[key[1]]
        mean = total / counts[key]
        matrix[i, j] = mean
        matrix[j, i] = mean
        if counts[key] > 1:
            duplicate_spreads.append(float(np.ptp(grouped_values[key])))

    unique_upper = len(positions) * (len(positions) + 1) // 2
    return matrix, {
        "source_values_in_crop": used,
        "unique_pair_cells": len(sums),
        "duplicate_pair_cells": int(sum(count > 1 for count in counts.values())),
        "maximum_duplicate_native_spread": max(duplicate_spreads, default=0.0),
        "upper_triangle_coverage_fraction": len(sums) / unique_upper,
    }


def reconstruct_poles(
    centers: np.ndarray, offsets: np.ndarray, predictions: np.ndarray, bin_bp: int
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int]]:
    centers = np.asarray(centers, dtype=np.int64)
    offsets = np.asarray(offsets, dtype=np.int64)
    positions = lattice_positions(centers, bin_bp)
    left = np.repeat(centers, len(offsets))
    right = (centers[:, None] + offsets[None, :]).reshape(-1)
    values = np.asarray(predictions, dtype=np.float64).reshape(-1)
    matrix, audit = aggregate_pairs(left, right, values, positions)
    return matrix, positions, audit


def reconstruct_chimaera(
    fragments: np.ndarray,
    starts: np.ndarray,
    distance_centers: np.ndarray,
    position_bin_bp: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int]]:
    # Chimaera's horizontal coordinate is the pair midpoint.  Its vertical
    # coordinate is genomic distance.  All released distance centers are even
    # multiples of the horizontal bin, so endpoints return to one common grid.
    left_all: list[np.ndarray] = []
    right_all: list[np.ndarray] = []
    value_all: list[np.ndarray] = []
    for fragment, start in zip(fragments, starts, strict=True):
        midpoints = int(start) + (np.arange(fragment.shape[1]) * position_bin_bp)
        midpoint_grid, distance_grid = np.meshgrid(
            midpoints, distance_centers, indexing="xy"
        )
        left_all.append((midpoint_grid - distance_grid // 2).reshape(-1))
        right_all.append((midpoint_grid + distance_grid // 2).reshape(-1))
        value_all.append(np.asarray(fragment, dtype=np.float64).reshape(-1))

    left = np.concatenate(left_all).astype(np.int64)
    right = np.concatenate(right_all).astype(np.int64)
    values = np.concatenate(value_all)
    positions = lattice_positions(np.concatenate([left, right]), position_bin_bp)
    matrix, audit = aggregate_pairs(left, right, values, positions)
    return matrix, positions, audit


def load_panels() -> tuple[list[SquarePanel], dict[str, dict[str, float | int]]]:
    panels: list[SquarePanel] = []
    audits: dict[str, dict[str, float | int]] = {}

    common = np.load(INPUTS, allow_pickle=False)
    observed_wt = np.asarray(common["observed_wt_crop_1kb"], dtype=np.float64)
    observed_deletion = np.asarray(common["observed_deletion_crop_1kb"], dtype=np.float64)
    panels.append(
        SquarePanel(
            "experimental", "Experimental Hi-TrAC", "1-kb pooled WT/deletion",
            "cLoops2 display: log10(PET + 1)", np.log10(observed_wt + 1),
            np.log10(observed_deletion + 1), EXPERIMENTAL_GRID_START,
            EXPERIMENTAL_GRID_END, 1_000,
            "Already a square 1-kb matrix; no geometric reconstruction.",
            "342 x 342 square raw-PET matrix", (PAPER_VMIN, PAPER_VMAX),
        )
    )

    with h5py.File(NATIVE / "akita_v2" / "preds.h5", "r") as handle:
        predictions = handle["preds"][:]
    wt, start, end = crop_square(square_from_upper(predictions[0, :, 0], 512, 2), 4_449_712, 2_048)
    deletion, _, _ = crop_square(square_from_upper(predictions[1, :, 0], 512, 2), 4_449_712, 2_048)
    panels.append(SquarePanel(
        "akita_v2", "AkitaV2", "HFF channel; 2,048-bp bins", "processed log(O/E)",
        wt, deletion, start, end, 2_048,
        "Released upper triangle mirrored into a square; masked diagonals remain missing.",
        "flattened upper triangle",
    ))

    z = np.load(NATIVE / "deepc" / "DEEPC_NATIVE_POLES.npz", allow_pickle=False)
    centers = np.rint(z["wt_coordinates"].mean(axis=1)).astype(np.int64)
    wt, positions, wt_audit = reconstruct_poles(centers, z["offsets_bp"], z["wt_predictions"], int(z["bin_bp"]))
    deletion, deletion_positions, deletion_audit = reconstruct_poles(
        centers, z["offsets_bp"], z["deletion_predictions"], int(z["bin_bp"])
    )
    if not np.array_equal(positions, deletion_positions):
        raise RuntimeError("DeepC WT/deletion coordinate lattices differ")
    audits["deepc_wt"] = wt_audit
    audits["deepc_deletion"] = deletion_audit
    panels.append(SquarePanel(
        "deepc", "DeepC", "K562 checkpoint; 5-kb bins", "normalized center-anchored interaction profile",
        wt, deletion, int(positions[0] - 2_500), int(positions[-1] + 2_500), 5_000,
        "Coordinate-equivalent implementation of Supplementary Figure 8: mirror poles, rotate to endpoint coordinates, crop, and average duplicate pair estimates; no value rescaling.",
        "83 center poles x 201 partner offsets",
    ))

    z = np.load(NATIVE / "orca" / "ORCA_NATIVE_MATRICES.npz", allow_pickle=False)
    wt, start, end = crop_square(z["hff_wt"], int(z["start"]), int(z["bin_bp"]))
    deletion, _, _ = crop_square(z["hff_deletion"], int(z["start"]), int(z["bin_bp"]))
    panels.append(SquarePanel(
        "orca", "Orca", "HFF channel; 4-kb bins", "normalized contact enrichment",
        wt, deletion, start, end, int(z["bin_bp"]), "Native square matrix cropped by genomic coordinate.", "dense symmetric matrix",
    ))

    z = np.load(NATIVE / "epcot" / "EPCOT_NATIVE_UPPER_TRIANGLE.npz", allow_pickle=False)
    n = int(z["matrix_bins"])
    wt, start, end = crop_square(square_from_upper(z["wt_control_dnase"], n, 0), int(z["start"]), int(z["bin_bp"]))
    deletion, _, _ = crop_square(square_from_upper(z["deletion_control_dnase"], n, 0), int(z["start"]), int(z["bin_bp"]))
    panels.append(SquarePanel(
        "epcot", "EPCOT", "HFF Micro-C head + K562 DNase; 1-kb bins", "model-native predicted O/E",
        wt, deletion, start, end, int(z["bin_bp"]), "Released upper triangle mirrored into a square.", "flattened upper triangle",
    ))

    z = np.load(NATIVE / "chromafold" / "CHROMAFOLD_NATIVE_VSTRIPES.npz", allow_pickle=False)
    centers = np.asarray(z["center_positions_derivative"], dtype=np.int64)
    wt, positions, wt_audit = reconstruct_poles(
        centers, z["partner_offsets_bp"], z["wt_control_dnase"], 10_000
    )
    deletion, deletion_positions, deletion_audit = reconstruct_poles(
        centers, z["partner_offsets_bp"], z["deletion_control_dnase"], 10_000
    )
    if not np.array_equal(positions, deletion_positions):
        raise RuntimeError("ChromaFold WT/deletion coordinate lattices differ")
    audits["chromafold_wt"] = wt_audit
    audits["chromafold_deletion"] = deletion_audit
    panels.append(SquarePanel(
        "chromafold", "ChromaFold motif", "K562 DNase proxy; 10-kb bins", "HiC-DC+ normalized Z-score",
        wt, deletion, int(positions[0]), int(positions[-1] + 10_000), 10_000,
        "Associated benchmark conversion: expand each V-stripe into genomic pairs and average duplicate (start, end) estimates; diagonal remains missing; no value rescaling.",
        "36 center V-stripes x 400 partner offsets",
    ))

    z = np.load(NATIVE / "alphagenome" / "ALPHAGENOME_NATIVE_TRACKDATA.npz", allow_pickle=False)
    wt, start, end = crop_square(z["reference_values"][:, :, 1], int(z["start"]), int(z["resolution"]))
    deletion, _, _ = crop_square(z["deletion_values"][:, :, 1], int(z["start"]), int(z["resolution"]))
    panels.append(SquarePanel(
        "alphagenome", "AlphaGenome", "HFFc6 Micro-C channel; 2,048-bp bins", "relative contact output",
        wt, deletion, start, end, int(z["resolution"]), "Native square matrix cropped by genomic coordinate.", "dense symmetric matrix",
    ))

    z = np.load(NATIVE / "chimaera" / "CHIMAERA_NATIVE_ROTATED_MAPS.npz", allow_pickle=False)
    wt, positions, wt_audit = reconstruct_chimaera(
        z["wt"], z["mapped_starts_derivative"], z["distance_centers_bp"], int(z["position_bin_bp"])
    )
    deletion, deletion_positions, deletion_audit = reconstruct_chimaera(
        z["deletion"], z["mapped_starts_derivative"], z["distance_centers_bp"], int(z["position_bin_bp"])
    )
    if not np.array_equal(positions, deletion_positions):
        raise RuntimeError("Chimaera WT/deletion coordinate lattices differ")
    audits["chimaera_wt"] = wt_audit
    audits["chimaera_deletion"] = deletion_audit
    panels.append(SquarePanel(
        "chimaera", "Chimaera", "generic human release; 2,048-bp endpoint grid", "standardized log-distance-residual",
        wt, deletion, int(positions[0] - 1_024), int(positions[-1] + 1_024), 2_048,
        "Invert native midpoint-distance coordinates to endpoint pairs and average overlapping fragments; unsupported distances remain missing; no interpolation or value rescaling.",
        "four overlapping 32 x 128 distance-position images",
    ))

    return panels, audits


def exact_limits(panel: SquarePanel) -> tuple[float, float, float]:
    if panel.fixed_limits is not None:
        low, high = panel.fixed_limits
    else:
        joined = np.concatenate(
            [panel.wt[np.isfinite(panel.wt)], panel.deletion[np.isfinite(panel.deletion)]]
        )
        low = float(joined.min())
        high = float(joined.max())
        if not high > low:
            high = low + 1.0
    difference = panel.deletion - panel.wt
    finite = np.abs(difference[np.isfinite(difference)])
    difference_max = float(finite.max()) if finite.size else 1.0
    return low, high, max(difference_max, np.finfo(float).eps)


def render(panels: list[SquarePanel]) -> None:
    configure_arial()
    red = mpl.colors.ListedColormap(PAPER_REDS, name="hitrac_paper_red")
    red.set_bad("#eeeeee")
    red.set_under("#ffffff")
    diverging = mpl.colormaps["RdBu_r"].copy()
    diverging.set_bad("#eeeeee")

    fig = plt.figure(figsize=(13.7, 2.42 * len(panels) + 1.25), constrained_layout=True)
    grid = fig.add_gridspec(
        len(panels), 5, width_ratios=[1, 1, 1, 0.045, 0.045], wspace=0.07, hspace=0.11
    )
    for row, panel in enumerate(panels):
        low, high, difference_max = exact_limits(panel)
        native_norm = mpl.colors.Normalize(low, high, clip=False)
        difference_norm = mpl.colors.TwoSlopeNorm(0, -difference_max, difference_max)
        arrays = [panel.wt, panel.deletion, panel.deletion - panel.wt]
        norms = [native_norm, native_norm, difference_norm]
        cmaps = [red, red, diverging]
        extent = [panel.start_bp / 1e6, panel.end_bp / 1e6, panel.end_bp / 1e6, panel.start_bp / 1e6]
        for column, (array, norm, cmap) in enumerate(zip(arrays, norms, cmaps, strict=True)):
            axis = fig.add_subplot(grid[row, column])
            axis.imshow(
                array, extent=extent, origin="upper", interpolation="nearest", aspect="equal",
                cmap=cmap, norm=norm, rasterized=True,
            )
            if row == 0:
                axis.set_title(["WT", "Inferred deletion", "Deletion - WT"][column], fontsize=10.5, fontweight="bold")
            if column == 0:
                axis.set_ylabel(
                    f"{panel.name}\n{panel.bin_bp:,}-bp bins", fontsize=8.2, fontweight="bold"
                )
            else:
                axis.set_yticklabels([])
            if row == len(panels) - 1:
                axis.set_xlabel("chr3 coordinate (Mb)", fontsize=8)
            else:
                axis.set_xticklabels([])
            axis.tick_params(labelsize=6.8, length=2, width=0.6)
            for spine in axis.spines.values():
                spine.set_linewidth(0.7)
            axis.text(
                0.01, 0.99, panel.scale, transform=axis.transAxes, va="top", ha="left",
                fontsize=6.1, bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none", "pad": 1.2},
            )

        native_bar_axis = fig.add_subplot(grid[row, 3])
        native_bar = fig.colorbar(
            mpl.cm.ScalarMappable(norm=native_norm, cmap=red), cax=native_bar_axis
        )
        native_bar.set_ticks([low, high])
        native_bar.ax.tick_params(labelsize=5.8, length=2, pad=1)
        native_bar.ax.set_title("native", fontsize=5.8, pad=2)

        difference_bar_axis = fig.add_subplot(grid[row, 4])
        difference_bar = fig.colorbar(
            mpl.cm.ScalarMappable(norm=difference_norm, cmap=diverging), cax=difference_bar_axis
        )
        difference_bar.set_ticks([-difference_max, difference_max])
        difference_bar.ax.tick_params(labelsize=5.8, length=2, pad=1)
        difference_bar.ax.set_title("delta", fontsize=5.8, pad=2)

    fig.suptitle(
        "BHLHE40 boundary-deletion benchmark: geometry-standardized square matrices\n"
        "chr3:4,803,502-5,144,387 (hg38); native values preserved; row scales are not comparable",
        fontsize=13, fontweight="bold",
    )
    for suffix in ("png", "svg"):
        fig.savefig(
            FIGURES / f"FIGURE3D_GEOMETRY_STANDARDIZED_SQUARES.{suffix}",
            dpi=320, bbox_inches="tight", facecolor="white",
        )
    plt.close(fig)


def write_outputs(
    panels: list[SquarePanel], reconstruction_audits: dict[str, dict[str, float | int]]
) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    audit_rows: list[dict[str, object]] = []
    for panel in panels:
        low, high, difference_max = exact_limits(panel)
        arrays[f"{panel.slug}_wt"] = panel.wt.astype(np.float32)
        arrays[f"{panel.slug}_deletion"] = panel.deletion.astype(np.float32)
        arrays[f"{panel.slug}_difference"] = (panel.deletion - panel.wt).astype(np.float32)
        arrays[f"{panel.slug}_start_bp"] = np.asarray(panel.start_bp, dtype=np.int64)
        arrays[f"{panel.slug}_end_bp"] = np.asarray(panel.end_bp, dtype=np.int64)
        arrays[f"{panel.slug}_bin_bp"] = np.asarray(panel.bin_bp, dtype=np.int64)
        finite = np.isfinite(panel.wt) & np.isfinite(panel.deletion)
        symmetry_error = max(
            float(np.nanmax(np.abs(panel.wt - panel.wt.T))),
            float(np.nanmax(np.abs(panel.deletion - panel.deletion.T))),
        )
        audit_rows.append(
            {
                "slug": panel.slug,
                "model": panel.name,
                "source_geometry": panel.source_geometry,
                "square_shape": list(panel.wt.shape),
                "bin_bp": panel.bin_bp,
                "start_bp": panel.start_bp,
                "end_bp": panel.end_bp,
                "native_scale": panel.scale,
                "display_min": low,
                "display_max": high,
                "difference_absolute_max": difference_max,
                "joint_finite_fraction": float(finite.mean()),
                "symmetry_max_abs_error": symmetry_error,
                "reconstruction": panel.reconstruction,
                "value_normalization_or_fitted_rescaling": False,
            }
        )

    common = np.load(INPUTS, allow_pickle=False)
    arrays["experimental_wt_raw_pet"] = np.asarray(common["observed_wt_crop_1kb"], dtype=np.float32)
    arrays["experimental_deletion_raw_pet"] = np.asarray(common["observed_deletion_crop_1kb"], dtype=np.float32)
    archive = OUT_DIR / "FIGURE3D_SQUARE_MATRICES.npz"
    np.savez_compressed(archive, **arrays)

    scale_table = OUT_DIR / "SQUARE_MATRIX_SCALE_AND_GEOMETRY.csv"
    with scale_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0].keys()))
        writer.writeheader()
        writer.writerows(audit_rows)

    audit = {
        "status": "PASS",
        "scientific_contract": {
            "standardized_property": "square endpoint-by-endpoint geometry over the same locus",
            "not_standardized": [
                "native bin size", "scientific target scale", "cell-type output head",
                "predicted distance range", "model-specific missing cells",
            ],
            "no_cross_model_value_normalization": True,
            "no_fitted_rescaling": True,
            "experimental_display": {
                "transform": "log10(PET + 1), matching the cLoops2/NAR-style display",
                "vmin": PAPER_VMIN,
                "vmax": PAPER_VMAX,
                "minimum_color": "white",
                "maximum_color": "red",
                "raw_pet_matrices_preserved_in_archive": True,
            },
        },
        "source_workflows": {
            "deepc": "C.Origami Supplementary Figure 8 geometry: mirror, rotate, crop; implemented equivalently by placing every pole value at its genomic endpoint pair.",
            "chromafold": "CBIGR/bulk_hic_benchmark bedpe_norm.py: expand 400-value V-stripes and average duplicate genomic pairs.",
            "chimaera": "Released midpoint-distance geometry inverted to endpoint pairs; overlapping fragment estimates averaged.",
        },
        "panels": audit_rows,
        "reconstruction_diagnostics": reconstruction_audits,
        "artifacts": {},
    }
    audit_path = OUT_DIR / "SQUARE_RECONSTRUCTION_AUDIT.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")

    artifacts = [
        archive,
        scale_table,
        audit_path,
        FIGURES / "FIGURE3D_GEOMETRY_STANDARDIZED_SQUARES.png",
        FIGURES / "FIGURE3D_GEOMETRY_STANDARDIZED_SQUARES.svg",
    ]
    hashes = [f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}" for path in artifacts]
    (OUT_DIR / "SQUARE_OUTPUT_SHA256SUMS.txt").write_text("\n".join(hashes) + "\n", encoding="utf-8")

    # Update the audit after all hash-stable primary artifacts exist.  The audit
    # intentionally does not include its own hash to avoid recursion.
    audit["artifacts"] = {
        path.relative_to(ROOT).as_posix(): sha256(path)
        for path in artifacts
        if path != audit_path
    }
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    hashes = [f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}" for path in artifacts]
    (OUT_DIR / "SQUARE_OUTPUT_SHA256SUMS.txt").write_text("\n".join(hashes) + "\n", encoding="utf-8")


def main() -> None:
    panels, reconstruction_audits = load_panels()
    render(panels)
    write_outputs(panels, reconstruction_audits)
    print(
        json.dumps(
            {
                "status": "PASS",
                "panels": [panel.name for panel in panels],
                "figure": str(FIGURES / "FIGURE3D_GEOMETRY_STANDARDIZED_SQUARES.png"),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
