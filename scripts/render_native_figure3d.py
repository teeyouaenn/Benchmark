#!/usr/bin/env python3
"""Render BHLHE40 WT/deletion perturbations in every model's native geometry."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "results" / "native"
INPUTS = ROOT / "results" / "inputs" / "COMMON_INPUTS_AND_TRUTH.npz"
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
DISPLAY_START = 4_803_501
DISPLAY_END = 5_144_387
DELETE_START = 4_976_067
DELETE_END = 4_976_790


@dataclass
class NativePanel:
    name: str
    subtitle: str
    wt: np.ndarray
    deletion: np.ndarray
    extent: tuple[float, float, float, float]
    xlabel: str = "chr3 coordinate (Mb)"
    ylabel: str = "chr3 coordinate (Mb)"
    aspect: str | float = "equal"


def square_from_upper(values: np.ndarray, n: int, diagonal: int) -> np.ndarray:
    matrix = np.full((n, n), np.nan, dtype=np.float32)
    indices = np.triu_indices(n, k=diagonal)
    if len(indices[0]) != len(values):
        raise ValueError(f"Cannot place {len(values)} values into n={n}, k={diagonal}")
    matrix[indices] = values
    matrix[(indices[1], indices[0])] = values
    return matrix


def crop_square(matrix: np.ndarray, start: int, bin_bp: int) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    centers = start + (np.arange(matrix.shape[0]) + 0.5) * bin_bp
    keep = (centers >= DISPLAY_START) & (centers < DISPLAY_END)
    indices = np.flatnonzero(keep)
    if len(indices) < 2:
        raise ValueError("Requested locus is outside model output")
    cropped = matrix[np.ix_(indices, indices)]
    lo = (start + indices[0] * bin_bp) / 1e6
    hi = (start + (indices[-1] + 1) * bin_bp) / 1e6
    return cropped, (lo, hi, hi, lo)


def robust_limits(wt: np.ndarray, deletion: np.ndarray) -> tuple[float, float, float]:
    joined = np.concatenate([wt[np.isfinite(wt)], deletion[np.isfinite(deletion)]])
    low = float(np.quantile(joined, 0.01))
    high = float(np.quantile(joined, 0.995))
    if not high > low:
        high = low + 1.0
    delta = deletion - wt
    delta_lim = float(np.quantile(np.abs(delta[np.isfinite(delta)]), 0.995))
    if delta_lim == 0:
        delta_lim = 1.0
    return low, high, delta_lim


def load_panels() -> tuple[NativePanel, list[NativePanel], list[NativePanel]]:
    common = np.load(INPUTS, allow_pickle=False)
    observed = NativePanel(
        "Experimental Hi-TrAC", "pooled WT/deletion, 1-kb log1p PET counts",
        np.log1p(common["observed_wt_crop_1kb"]), np.log1p(common["observed_deletion_crop_1kb"]),
        (DISPLAY_START / 1e6, DISPLAY_END / 1e6, DISPLAY_END / 1e6, DISPLAY_START / 1e6),
    )

    dense: list[NativePanel] = []
    with h5py.File(NATIVE / "akita_v2" / "preds.h5", "r") as handle:
        pred = handle["preds"][:]
    akita_wt = square_from_upper(pred[0, :, 0], 512, 2)
    akita_del = square_from_upper(pred[1, :, 0], 512, 2)
    akita_wt, extent = crop_square(akita_wt, 4_449_712, 2048)
    akita_del, _ = crop_square(akita_del, 4_449_712, 2048)
    dense.append(NativePanel("AkitaV2", "HFF channel; 2,048-bp processed log(O/E)", akita_wt, akita_del, extent))

    z = np.load(NATIVE / "orca" / "ORCA_NATIVE_MATRICES.npz", allow_pickle=False)
    wt, extent = crop_square(z["hff_wt"], int(z["start"]), int(z["bin_bp"]))
    deletion, _ = crop_square(z["hff_deletion"], int(z["start"]), int(z["bin_bp"]))
    dense.append(NativePanel("Orca", "HFF channel; 4-kb normalized contact enrichment", wt, deletion, extent))

    z = np.load(NATIVE / "epcot" / "EPCOT_NATIVE_UPPER_TRIANGLE.npz", allow_pickle=False)
    epcot_wt = square_from_upper(z["wt_control_dnase"], int(z["matrix_bins"]), 0)
    epcot_del = square_from_upper(z["deletion_control_dnase"], int(z["matrix_bins"]), 0)
    wt, extent = crop_square(epcot_wt, int(z["start"]), int(z["bin_bp"]))
    deletion, _ = crop_square(epcot_del, int(z["start"]), int(z["bin_bp"]))
    dense.append(NativePanel("EPCOT", "HFF Micro-C head + K562 control DNase; 1-kb native scale", wt, deletion, extent))

    z = np.load(NATIVE / "alphagenome" / "ALPHAGENOME_NATIVE_TRACKDATA.npz", allow_pickle=False)
    wt, extent = crop_square(z["reference_values"][:, :, 1], int(z["start"]), int(z["resolution"]))
    deletion, _ = crop_square(z["deletion_values"][:, :, 1], int(z["start"]), int(z["resolution"]))
    dense.append(NativePanel("AlphaGenome", "HFFc6 Micro-C channel; 2,048-bp relative contact output", wt, deletion, extent))

    nonsquare: list[NativePanel] = []
    z = np.load(NATIVE / "deepc" / "DEEPC_NATIVE_POLES.npz", allow_pickle=False)
    centers = z["wt_coordinates"].mean(axis=1)
    keep = (centers >= DISPLAY_START) & (centers < DISPLAY_END)
    wt = z["wt_predictions"][keep]
    deletion = z["deletion_predictions"][keep]
    nonsquare.append(NativePanel(
        "DeepC", "K562 5-kb native center poles",
        wt, deletion,
        (z["offsets_bp"][0] / 1000, z["offsets_bp"][-1] / 1000,
         centers[keep][-1] / 1e6, centers[keep][0] / 1e6),
        "partner offset from center (kb)", "center coordinate (Mb)", "auto",
    ))

    z = np.load(NATIVE / "chromafold" / "CHROMAFOLD_NATIVE_VSTRIPES.npz", allow_pickle=False)
    centers = z["center_positions_reference_wt"]
    keep = (centers >= DISPLAY_START) & (centers < DISPLAY_END)
    nonsquare.append(NativePanel(
        "ChromaFold motif", "10-kb V-stripes; K562 DNase as scATAC proxy",
        z["wt_control_dnase"][keep], z["deletion_control_dnase"][keep],
        (z["partner_offsets_bp"][0] / 1000, z["partner_offsets_bp"][-1] / 1000,
         centers[keep][-1] / 1e6, centers[keep][0] / 1e6),
        "partner offset from center (kb)", "center coordinate (Mb)", "auto",
    ))

    z = np.load(NATIVE / "chimaera" / "CHIMAERA_NATIVE_ROTATED_MAPS.npz", allow_pickle=False)
    wt = np.concatenate(list(z["wt"]), axis=1)
    deletion = np.concatenate(list(z["deletion"]), axis=1)
    nonsquare.append(NativePanel(
        "Chimaera", "four overlapping native rotated images; 2,048-bp position axis",
        wt, deletion,
        (0, wt.shape[1], z["distance_centers_bp"][-1] / 1000, z["distance_centers_bp"][0] / 1000),
        "concatenated local position bins", "genomic separation (kb)", "auto",
    ))
    return observed, dense, nonsquare


def decorate_axis(ax: mpl.axes.Axes, panel: NativePanel, column: int) -> None:
    ax.set_xlabel(panel.xlabel, fontsize=7)
    if column == 0:
        ax.set_ylabel(f"{panel.name}\n{panel.ylabel}", fontsize=8, fontweight="bold")
    else:
        ax.set_ylabel("")
    ax.tick_params(labelsize=6, length=2)


def draw_panels(panels: list[NativePanel], path: Path, title: str) -> None:
    fig, axes = plt.subplots(len(panels), 3, figsize=(11.5, 3.0 * len(panels)), constrained_layout=True)
    if len(panels) == 1:
        axes = np.asarray([axes])
    for row, panel in enumerate(panels):
        low, high, dlim = robust_limits(panel.wt, panel.deletion)
        arrays = [panel.wt, panel.deletion, panel.deletion - panel.wt]
        cmaps = ["Reds", "Reds", "RdBu_r"]
        norms = [mpl.colors.Normalize(low, high), mpl.colors.Normalize(low, high), mpl.colors.TwoSlopeNorm(0, -dlim, dlim)]
        labels = ["WT", "inferred deletion", "deletion - WT"]
        for col, (array, cmap, norm, label) in enumerate(zip(arrays, cmaps, norms, labels)):
            ax = axes[row, col]
            image = ax.imshow(array, cmap=cmap, norm=norm, extent=panel.extent, interpolation="nearest", aspect=panel.aspect)
            if row == 0:
                ax.set_title(label, fontsize=11, fontweight="bold")
            decorate_axis(ax, panel, col)
            if col == 2:
                cb = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
                cb.ax.tick_params(labelsize=6)
            ax.text(0.01, 0.99, panel.subtitle, transform=ax.transAxes, va="top", ha="left",
                    fontsize=6.5, color="black", bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none", "pad": 1.5})
    fig.suptitle(title + "\nchr3:4,803,502-5,144,387 (hg38); inferred 723-bp BHLHE40-boundary deletion",
                 fontsize=13, fontweight="bold")
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def draw_all(panels: list[NativePanel], path: Path) -> None:
    fig, axes = plt.subplots(len(panels), 3, figsize=(12, 2.45 * len(panels)), constrained_layout=True)
    for row, panel in enumerate(panels):
        low, high, dlim = robust_limits(panel.wt, panel.deletion)
        for col, (array, cmap, norm) in enumerate([
            (panel.wt, "Reds", mpl.colors.Normalize(low, high)),
            (panel.deletion, "Reds", mpl.colors.Normalize(low, high)),
            (panel.deletion - panel.wt, "RdBu_r", mpl.colors.TwoSlopeNorm(0, -dlim, dlim)),
        ]):
            ax = axes[row, col]
            ax.imshow(array, cmap=cmap, norm=norm, extent=panel.extent, interpolation="nearest", aspect=panel.aspect)
            if row == 0:
                ax.set_title(["WT", "inferred deletion", "deletion - WT"][col], fontsize=11, fontweight="bold")
            decorate_axis(ax, panel, col)
            ax.text(0.01, 0.99, panel.subtitle, transform=ax.transAxes, va="top", fontsize=6,
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1})
    fig.suptitle("BHLHE40 boundary perturbation in seven predictors' native outputs\n"
                 "Rows retain native resolution, geometry, channel and scientific scale; colors are row-specific",
                 fontsize=13, fontweight="bold")
    fig.savefig(path, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def draw_target_assisted(path: Path) -> None:
    rows = []
    z = np.load(NATIVE / "epcot" / "EPCOT_NATIVE_UPPER_TRIANGLE.npz", allow_pickle=False)
    n = int(z["matrix_bins"])
    matrices = [square_from_upper(z[k], n, 0) for k in
                ["wt_control_dnase", "deletion_control_dnase", "deletion_target_assisted_hitrac1d"]]
    cropped = [crop_square(x, int(z["start"]), int(z["bin_bp"]))[0] for x in matrices]
    extent = crop_square(matrices[0], int(z["start"]), int(z["bin_bp"]))[1]
    rows.append(("EPCOT", "1-kb HFF Micro-C head", cropped, extent, "equal"))

    z = np.load(NATIVE / "chromafold" / "CHROMAFOLD_NATIVE_VSTRIPES.npz", allow_pickle=False)
    centers = z["center_positions_reference_wt"]
    keep = (centers >= DISPLAY_START) & (centers < DISPLAY_END)
    arrays = [z[k][keep] for k in ["wt_control_dnase", "deletion_control_dnase", "deletion_hitrac1d_assisted"]]
    extent = (z["partner_offsets_bp"][0] / 1000, z["partner_offsets_bp"][-1] / 1000,
              centers[keep][-1] / 1e6, centers[keep][0] / 1e6)
    rows.append(("ChromaFold motif", "10-kb V-stripes", arrays, extent, "auto"))

    fig, axes = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True)
    for row, (name, subtitle, arrays, extent, aspect) in enumerate(rows):
        joined = np.concatenate([a[np.isfinite(a)] for a in arrays])
        low, high = np.quantile(joined, [0.01, 0.995])
        delta = arrays[2] - arrays[1]
        dlim = max(float(np.quantile(np.abs(delta[np.isfinite(delta)]), 0.995)), 1e-6)
        for col, array in enumerate(arrays + [delta]):
            norm = mpl.colors.Normalize(low, high) if col < 3 else mpl.colors.TwoSlopeNorm(0, -dlim, dlim)
            axes[row, col].imshow(array, cmap="Reds" if col < 3 else "RdBu_r", norm=norm,
                                  extent=extent, interpolation="nearest", aspect=aspect)
            if row == 0:
                axes[row, col].set_title(["WT + control DNase", "deletion + shifted control DNase",
                                          "deletion + Hi-TrAC 1D", "assisted - strict deletion"][col], fontsize=9)
            if col == 0:
                axes[row, col].set_ylabel(f"{name}\n{subtitle}", fontweight="bold", fontsize=8)
            axes[row, col].tick_params(labelsize=6)
    fig.suptitle("Target-assisted sensitivity analysis — not de novo and ineligible for model ranking\n"
                 "Experimental deletion Hi-TrAC 1D was quantile-matched to control-DNase units",
                 fontsize=12, fontweight="bold", color="#8b0000")
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_perturbation_summary(panels: list[NativePanel]) -> None:
    rows = []
    for panel in panels:
        mask = np.isfinite(panel.wt) & np.isfinite(panel.deletion)
        wt = panel.wt[mask].astype(float)
        deletion = panel.deletion[mask].astype(float)
        delta = deletion - wt
        pearson = float(np.corrcoef(wt, deletion)[0, 1]) if wt.size > 1 else float("nan")
        rows.append({
            "model": panel.name,
            "native_values_compared": int(wt.size),
            "wt_deletion_pearson": pearson,
            "mean_absolute_native_delta": float(np.mean(np.abs(delta))),
            "rms_native_delta": float(np.sqrt(np.mean(delta**2))),
            "rms_delta_over_wt_sd": float(np.sqrt(np.mean(delta**2)) / (np.std(wt) + 1e-12)),
            "warning": "Within-model perturbation magnitude only; native scales are not comparable across rows.",
        })
    (RESULTS / "NATIVE_PERTURBATION_SUMMARY.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")


def draw_locus_and_inputs(path: Path) -> None:
    common = np.load(INPUTS, allow_pickle=False)
    genes = json.loads((ROOT / "manifests" / "gene_track.json").read_text(encoding="utf-8"))["genes"]
    lo = DISPLAY_START / 1e6
    hi = DISPLAY_END / 1e6
    x = (3_000_000 + (np.arange(4000) + 0.5) * 1000) / 1e6
    keep = (x >= lo) & (x < hi)
    fig, axes = plt.subplots(4, 1, figsize=(12, 6.8), sharex=True,
                             gridspec_kw={"height_ratios": [0.9, 1, 1, 1]})
    gene_ax = axes[0]
    gene_ax.set_ylim(-0.2, 1.8)
    for lane, gene in enumerate(genes):
        start = max(gene["start"] / 1e6, lo)
        end = min(gene["end"] / 1e6, hi)
        if end <= start:
            continue
        y = 0.35 if lane % 2 == 0 else 1.15
        gene_ax.plot([start, end], [y, y], lw=4, color="#3b82f6")
        gene_ax.text((start + end) / 2, y + 0.18, gene["name"], ha="center", va="bottom",
                     fontsize=9, fontstyle="italic", color="#1d4ed8")
        arrow_x = end if gene["strand"] > 0 else start
        gene_ax.scatter([arrow_x], [y], marker=">" if gene["strand"] > 0 else "<", s=45, color="#1d4ed8")
    gene_ax.axvspan(DELETE_START / 1e6, DELETE_END / 1e6, color="black", alpha=0.85)
    gene_ax.text((DELETE_START + DELETE_END) / 2 / 1e6, 1.65, "inferred 723-bp deletion", ha="center", fontsize=8)
    gene_ax.axis("off")

    tracks = [
        (common["control_dnase_1kb"], "WT input: control K562 DNase", "#7e22ce", "DNase"),
        (common["deletion_control_dnase_1kb"], "Deletion input: same control DNase shifted with derivative DNA", "#7e22ce", "DNase"),
        (common["target_assisted_surrogate_1kb_derivative"],
         "Target-assisted input: deletion Hi-TrAC 1D rank-matched to DNase — NOT DE NOVO", "#b91c1c", "Hi-TrAC 1D"),
    ]
    ymax = max(float(np.quantile(values[keep], 0.995)) for values, _, _, _ in tracks)
    for ax, (values, label, color, ylabel) in zip(axes[1:], tracks):
        ax.fill_between(x[keep], values[keep], color=color, alpha=0.85, linewidth=0)
        ax.set_ylim(0, ymax * 1.05)
        ax.set_ylabel(f"{ylabel}\n(shared units)", fontsize=8)
        ax.set_title(label, loc="left", fontsize=9)
        ax.axvspan(DELETE_START / 1e6, DELETE_END / 1e6, color="black", alpha=0.25)
        ax.grid(axis="y", alpha=0.15)
    axes[-1].set_xlim(lo, hi)
    axes[-1].set_xlabel("chr3 coordinate (Mb)")
    fig.suptitle("BHLHE40 Figure 3D locus and accessibility inputs", fontsize=13, fontweight="bold", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    observed, dense, nonsquare = load_panels()
    draw_panels([observed] + dense, FIGURES / "FIGURE3D_DENSE_NATIVE_OUTPUTS.png",
                "NAR Figure 3D locus — experimental map and native square-map predictors")
    draw_panels(nonsquare, FIGURES / "FIGURE3D_NON_SQUARE_NATIVE_OUTPUTS.png",
                "NAR Figure 3D locus — native pole, V-stripe and rotated-image predictors")
    draw_all([observed] + dense + nonsquare, FIGURES / "FIGURE3D_ALL_NATIVE_OUTPUTS.png")
    draw_target_assisted(FIGURES / "FIGURE3D_TARGET_ASSISTED_NOT_DENOVO.png")
    draw_locus_and_inputs(FIGURES / "FIGURE3D_LOCUS_AND_INPUT_TRACKS.png")
    write_perturbation_summary([observed] + dense + nonsquare)
    print(json.dumps({"status": "COMPLETE", "figures": sorted(str(p) for p in FIGURES.glob("*.png"))}, indent=2))


if __name__ == "__main__":
    main()
