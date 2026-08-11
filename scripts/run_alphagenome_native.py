#!/usr/bin/env python3
"""Run official AlphaGenome all-fold weights on the BHLHE40 deletion."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import jax
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome_research.model import dna_model


DELETE_START = 4_976_067  # 0-based, inclusive
DELETE_END = 4_976_790  # 0-based, exclusive
INTERVAL_START = 4_449_712
INTERVAL_END = 5_498_288


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--common-npz", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--gpu", default="0")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    devices = jax.devices("gpu")
    if not devices:
        raise RuntimeError("AlphaGenome requires a visible GPU")
    device = devices[int(args.gpu)] if len(devices) > int(args.gpu) else devices[0]
    model = dna_model.create_from_huggingface("all_folds", device=device)

    metadata = dna_model.metadata_lib.load(dna_model.Organism.HOMO_SAPIENS)
    contact_metadata = metadata.contact_maps.copy()
    if contact_metadata is None or contact_metadata.empty:
        raise RuntimeError("released model contains no human contact-map metadata")
    ontology_terms = sorted(
        term
        for term in contact_metadata["ontology_curie"].dropna().unique().tolist()
        if str(term).strip()
    )

    common = np.load(args.common_npz)
    wt = np.asarray(common["wt_sequence"], dtype=np.uint8).tobytes().decode("ascii")
    common_start = int(common["start"])
    rel_start = DELETE_START - common_start
    rel_end = DELETE_END - common_start
    anchor_base = wt[rel_start - 1]
    deleted_bases = wt[rel_start:rel_end]
    reference_allele = anchor_base + deleted_bases
    alternate_allele = anchor_base

    # AlphaGenome's Variant follows VCF-style 1-based positioning.  The anchor
    # is the base immediately before the 0-based half-open deleted interval.
    variant = genome.Variant(
        chromosome="chr3",
        position=DELETE_START,
        reference_bases=reference_allele,
        alternate_bases=alternate_allele,
    )
    interval = genome.Interval(
        chromosome="chr3", start=INTERVAL_START, end=INTERVAL_END
    )
    output = model.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=[dna_model.OutputType.CONTACT_MAPS],
        ontology_terms=ontology_terms,
    )
    reference = output.reference.contact_maps
    alternate = output.alternate.contact_maps
    ref_values = np.asarray(reference.values, dtype=np.float32)
    alt_values = np.asarray(alternate.values, dtype=np.float32)
    if ref_values.shape != alt_values.shape or ref_values.ndim != 3:
        raise RuntimeError(
            f"unexpected AlphaGenome shapes: {ref_values.shape}, {alt_values.shape}"
        )
    if ref_values.shape[:2] != (512, 512):
        raise RuntimeError(f"unexpected contact-map geometry: {ref_values.shape}")
    if ref_values.shape[2] == 0:
        raise RuntimeError("ontology filtering returned zero contact-map tracks")
    if not np.isfinite(ref_values).all() or not np.isfinite(alt_values).all():
        raise RuntimeError("non-finite AlphaGenome contact prediction")

    native_npz = args.out_dir / "ALPHAGENOME_NATIVE_TRACKDATA.npz"
    np.savez_compressed(
        native_npz,
        reference_values=ref_values,
        deletion_values=alt_values,
        start=np.int64(reference.interval.start),
        end=np.int64(reference.interval.end),
        resolution=np.int64(reference.resolution),
    )
    metadata_csv = args.out_dir / "ALPHAGENOME_CONTACT_METADATA.csv"
    reference.metadata.to_csv(metadata_csv, index=False)
    audit_path = args.out_dir / "ALPHAGENOME_RUN_AUDIT.json"
    audit = {
        "status": "COMPLETE",
        "model": "official AlphaGenome all-fold local weights",
        "input": "DNA only",
        "k562_matched": False,
        "caveat": "released human contact-map metadata contains no K562 track",
        "variant": {
            "chromosome": "chr3",
            "position_1based_anchor": DELETE_START,
            "deleted_interval_0based_half_open": [DELETE_START, DELETE_END],
            "reference_allele_length": len(reference_allele),
            "alternate_allele_length": len(alternate_allele),
        },
        "requested_interval_0based_half_open": [INTERVAL_START, INTERVAL_END],
        "returned_interval_0based_half_open": [
            int(reference.interval.start),
            int(reference.interval.end),
        ],
        "native_resolution_bp": int(reference.resolution),
        "native_shape": list(ref_values.shape),
        "ontology_terms_requested": ontology_terms,
        "tracks_returned": int(ref_values.shape[2]),
        "common_input_sha256": sha256(args.common_npz),
        "native_npz_sha256": sha256(native_npz),
        "metadata_csv_sha256": sha256(metadata_csv),
        "summary": {
            "reference": {
                "minimum": float(ref_values.min()),
                "maximum": float(ref_values.max()),
                "mean": float(ref_values.mean()),
                "standard_deviation": float(ref_values.std()),
            },
            "deletion": {
                "minimum": float(alt_values.min()),
                "maximum": float(alt_values.max()),
                "mean": float(alt_values.mean()),
                "standard_deviation": float(alt_values.std()),
            },
        },
    }
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
