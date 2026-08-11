#!/usr/bin/env python3
"""Render the historical K562.hic locus and combine it with Gviz tracks."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.image import imread


BIN_BP = 1000
MATRIX_START = 4_803_000
MATRIX_END = 5_145_000
LOCUS_START = 4_803_502
LOCUS_END = 5_144_387
DELETION_START = 4_976_067
DELETION_END = 4_976_790
SUBTAD_START = 4_977_047
SUBTAD_END = 5_056_047


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_dump(path: Path) -> tuple[np.ndarray, list[tuple[int, int, float]]]:
    bins = np.arange(MATRIX_START, MATRIX_END, BIN_BP, dtype=np.int64)
    index = {int(position): i for i, position in enumerate(bins)}
    matrix = np.zeros((bins.size, bins.size), dtype=np.float64)
    records: list[tuple[int, int, float]] = []
    for line in path.read_text(encoding="ascii", errors="ignore").splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            left, right = int(fields[0]), int(fields[1])
            value = float(fields[2])
        except ValueError:
            continue
        if left not in index or right not in index:
            continue
        i, j = index[left], index[right]
        matrix[i, j] += value
        if i != j:
            matrix[j, i] += value
        records.append((left, right, value))
    if not records:
        raise RuntimeError("no valid contact records were parsed")
    return matrix, records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--ideogram", type=Path, required=True)
    parser.add_argument("--locus-tracks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--clean-tsv", type=Path, required=True)
    parser.add_argument("--audit-json", type=Path, required=True)
    args = parser.parse_args()

    matrix, records = parse_dump(args.dump)
    args.clean_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.clean_tsv.write_text(
        "bin1_start_0based\tbin2_start_0based\tobserved_count\n"
        + "".join(f"{a}\t{b}\t{v:.12g}\n" for a, b, v in records),
        encoding="ascii",
    )

    upper_mass = float(np.triu(matrix).sum())
    diagonal_mass = float(np.diag(matrix).sum())
    nonzero_upper = int(np.count_nonzero(np.triu(matrix)))
    row_marginal = matrix.sum(axis=1)
    transformed = np.log10(1.0 + matrix)
    extent = [MATRIX_START / 1e6, MATRIX_END / 1e6,
              MATRIX_END / 1e6, MATRIX_START / 1e6]
    red_white = LinearSegmentedColormap.from_list(
        "paper_red", ["#ffffff", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"]
    )

    ideogram = imread(args.ideogram)
    locus_tracks = imread(args.locus_tracks)
    figure = plt.figure(figsize=(12, 14), constrained_layout=False)
    grid = figure.add_gridspec(
        4, 1, height_ratios=[0.9, 1.75, 1.0, 6.4], hspace=0.18,
        left=0.085, right=0.90, top=0.975, bottom=0.06
    )

    axis = figure.add_subplot(grid[0])
    axis.imshow(ideogram)
    axis.axis("off")
    axis.text(-0.035, 0.85, "A", transform=axis.transAxes, fontsize=18, fontweight="bold")

    axis = figure.add_subplot(grid[1])
    axis.imshow(locus_tracks)
    axis.axis("off")
    axis.text(-0.035, 0.95, "B", transform=axis.transAxes, fontsize=18, fontweight="bold")

    x_mb = (np.arange(matrix.shape[0]) * BIN_BP + MATRIX_START + BIN_BP / 2) / 1e6
    axis = figure.add_subplot(grid[2])
    axis.fill_between(x_mb, row_marginal, color="#7a2cb8", alpha=0.95, linewidth=0)
    axis.plot(x_mb, row_marginal, color="#6a1b9a", linewidth=0.7)
    axis.set_xlim(MATRIX_START / 1e6, MATRIX_END / 1e6)
    axis.set_ylabel("1D contact\nmarginal", fontsize=10)
    axis.set_xticklabels([])
    axis.spines[["top", "right"]].set_visible(False)
    axis.axvspan(DELETION_START / 1e6, DELETION_END / 1e6, color="#d00000", alpha=0.35)
    axis.text(-0.035, 0.88, "C", transform=axis.transAxes, fontsize=18, fontweight="bold")

    axis = figure.add_subplot(grid[3])
    image = axis.imshow(
        transformed,
        origin="upper",
        extent=extent,
        interpolation="nearest",
        cmap=red_white,
        vmin=0.0,
        vmax=2.0,
        aspect="equal",
        rasterized=True,
    )
    axis.set_xlabel("chr3 coordinate (Mb)", fontsize=11)
    axis.set_ylabel("chr3 coordinate (Mb)", fontsize=11)
    axis.set_title(
        "Historical K562.hic - observed NONE, native 1-kb bins\n"
        "white-red scale = log10(observed contact + 1), fixed 0-2",
        fontsize=12,
        fontweight="bold",
    )
    axis.plot(
        [SUBTAD_START / 1e6, SUBTAD_END / 1e6, SUBTAD_END / 1e6,
         SUBTAD_START / 1e6, SUBTAD_START / 1e6],
        [SUBTAD_START / 1e6, SUBTAD_START / 1e6, SUBTAD_END / 1e6,
         SUBTAD_END / 1e6, SUBTAD_START / 1e6],
        linestyle="--", color="#1565c0", linewidth=1.7,
        label="published K562 active sub-TAD"
    )
    deletion_mid = (DELETION_START + DELETION_END) / 2 / 1e6
    axis.axvline(deletion_mid, color="#d00000", linewidth=1.1, alpha=0.9)
    axis.axhline(deletion_mid, color="#d00000", linewidth=1.1, alpha=0.9)
    axis.legend(loc="lower right", frameon=True, fontsize=8)
    axis.text(-0.035, 0.98, "D", transform=axis.transAxes, fontsize=18, fontweight="bold")
    color_axis = figure.add_axes([0.92, 0.075, 0.022, 0.37])
    colorbar = figure.colorbar(image, cax=color_axis)
    colorbar.set_label("log10(contact + 1)", fontsize=10)

    figure.suptitle(
        "Audited BHLHE40 Figure 3D locus and historical K562 contact map",
        fontsize=15,
        fontweight="bold",
        y=0.997,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=240, facecolor="white")
    plt.close(figure)

    audit = {
        "hic_source": {
            "path": "/data/vinhtb/activhitracMIX_stage/data/user_k562_hitrac_mix/K562.hic",
            "sha256": "8ae9a61bf8f420d33f05f8c32747e520be57087650f1e6d0b499c2893673f9b3",
            "extraction": "Juicer dump observed NONE, BP 1000",
        },
        "requested_locus_1based_inclusive": [LOCUS_START, LOCUS_END],
        "enclosing_hic_grid_0based_half_open": [MATRIX_START, MATRIX_END],
        "matrix_shape": list(matrix.shape),
        "upper_triangle_contact_mass": upper_mass,
        "diagonal_contact_mass": diagonal_mass,
        "upper_triangle_nonzero_cells": nonzero_upper,
        "display_transform": "log10(observed_count + 1), fixed range 0 to 2",
        "clean_dump_sha256": sha256(args.clean_tsv),
        "figure_sha256": sha256(args.output),
        "lineage_warning": (
            "This historical K562.hic is not the pooled/downsampled experimental "
            "control used in NAR Figure 3D. It is shown because the user requested "
            "their K562.hic file for the identical genomic locus."
        ),
    }
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)
    args.audit_json.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
