#!/usr/bin/env python3
"""Run the released human Chimaera DNA encoder and Hi-C decoder.

The output remains in Chimaera's native rotated distance-coordinate geometry
(32 distance rows by 128 genomic-position columns) and native standardized
log-distance-residual scale.  Forward and reverse-complement predictions are
averaged exactly as in the released ``ModelContainer.predict`` method.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import torch


CHROM = "chr3"
COMMON_START = 3_000_000
DISPLAY_START = 4_803_501
DISPLAY_END = 5_144_387
DELETE_START = 4_976_067
DELETE_END = 4_976_790
DELETE_LEN = DELETE_END - DELETE_START
DNA_BP = 524_288
MAPPED_BP = 262_144
OFFSET_BP = 131_072
POSITION_BIN_BP = 2_048
DISTANCE_STEP_BP = 4_096


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


def one_hot(sequence: np.ndarray) -> torch.Tensor:
    """One-hot encode ASCII A/C/G/T in Chimaera's released channel order."""
    sequence = np.char.upper(sequence.astype("U1"))
    encoded = np.zeros((4, sequence.size), dtype=np.float32)
    for channel, base in enumerate(("A", "C", "G", "T")):
        encoded[channel] = sequence == base
    return torch.from_numpy(encoded)


def rc_tensor(x: torch.Tensor) -> torch.Tensor:
    # A,C,G,T -> T,G,C,A and reverse genomic order.
    return torch.flip(x, dims=[1, 2])


def infer(
    dna_encoder: torch.jit.ScriptModule,
    hic_decoder: torch.jit.ScriptModule,
    x: torch.Tensor,
    device: torch.device,
) -> np.ndarray:
    with torch.inference_mode():
        x = x.to(device)
        forward = hic_decoder(dna_encoder(x)[0])
        reverse = hic_decoder(dna_encoder(rc_tensor(x))[0])
        # Released code flips the genomic-position dimension of the RC output.
        reverse = torch.flip(reverse, dims=[3])
        output = 0.5 * (forward + reverse)
    output = output[:, 0].cpu().numpy().astype(np.float32)
    if output.shape[1:] != (32, 128):
        raise RuntimeError(f"unexpected Chimaera output shape {output.shape}")
    if not np.isfinite(output).all():
        raise RuntimeError("non-finite Chimaera prediction")
    return output


def extract_windows(sequence: np.ndarray, mapped_starts: np.ndarray) -> torch.Tensor:
    windows = []
    for mapped_start in mapped_starts:
        dna_start = int(mapped_start) - OFFSET_BP
        begin = dna_start - COMMON_START
        end = begin + DNA_BP
        if begin < 0 or end > sequence.size:
            raise RuntimeError("Chimaera DNA window exceeds registered common sequence")
        windows.append(one_hot(sequence[begin:end]))
    return torch.stack(windows, dim=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--common-inputs", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--device", default="cuda:1")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    dna_encoder_path = args.model_dir / "dna_encoder.pt"
    hic_decoder_path = args.model_dir / "hic_decoder.pt"
    dna_encoder = torch.jit.load(str(dna_encoder_path), map_location=device).eval()
    hic_decoder = torch.jit.load(str(hic_decoder_path), map_location=device).eval()

    common = np.load(args.common_inputs)
    wt_sequence = np.asarray(common["wt_sequence"], dtype=np.uint8).view("S1")
    deletion_sequence = np.asarray(common["deletion_sequence"], dtype=np.uint8).view("S1")

    first = (DISPLAY_START // OFFSET_BP) * OFFSET_BP
    mapped_starts = np.arange(first, DISPLAY_END, OFFSET_BP, dtype=np.int64)
    if mapped_starts[-1] + MAPPED_BP < DISPLAY_END:
        mapped_starts = np.append(mapped_starts, mapped_starts[-1] + OFFSET_BP)

    wt_x = extract_windows(wt_sequence, mapped_starts)
    deletion_x = extract_windows(deletion_sequence, mapped_starts)
    wt = infer(dna_encoder, hic_decoder, wt_x, device)
    deletion = infer(dna_encoder, hic_decoder, deletion_x, device)

    reference_mapped_starts_deletion = np.where(
        mapped_starts < DELETE_START, mapped_starts, mapped_starts + DELETE_LEN
    )
    distance_centers_bp = (
        np.arange(32, dtype=np.int64) + 2
    ) * DISTANCE_STEP_BP
    output = args.out_dir / "CHIMAERA_NATIVE_ROTATED_MAPS.npz"
    np.savez_compressed(
        output,
        wt=wt,
        deletion=deletion,
        difference=deletion - wt,
        mapped_starts_derivative=mapped_starts,
        mapped_starts_reference_wt=mapped_starts,
        mapped_starts_reference_deletion=reference_mapped_starts_deletion,
        mapped_span_bp=np.int64(MAPPED_BP),
        position_bin_bp=np.int64(POSITION_BIN_BP),
        distance_centers_bp=distance_centers_bp,
    )

    audit = {
        "status": "COMPLETE",
        "scientific_status": "OFFICIAL_HUMAN_RELEASE_NO_K562_SPECIFIC_CHANNEL",
        "source": {"git_head": git_head(args.source)},
        "checkpoint": {
            "dna_encoder_sha256": sha256(dna_encoder_path),
            "hic_decoder_sha256": sha256(hic_decoder_path),
            "model_params_sha256": sha256(args.model_dir / "model_params.json"),
            "data_params_sha256": sha256(args.model_dir / "data_params.json"),
        },
        "native_prediction": {
            "description": "standardized log-distance-residual contact signal",
            "geometry": "rotated distance-coordinate image",
            "shape_per_condition": list(wt.shape),
            "fragment_shape": [32, 128],
            "dna_input_bp": DNA_BP,
            "mapped_span_bp": MAPPED_BP,
            "position_bin_bp": POSITION_BIN_BP,
            "distance_step_bp": DISTANCE_STEP_BP,
            "output_file": str(output),
            "output_sha256": sha256(output),
        },
        "inference": {
            "forward_reverse_complement_average": True,
            "implementation_equivalence": (
                "direct TorchScript composition of released dna_encoder.pt and "
                "hic_decoder.pt; identical composition and RC averaging to the "
                "released ModelContainer"
            ),
        },
        "deletion_status": (
            "inferred clean SpCas9 cut-to-cut deletion; exact clone junction unavailable"
        ),
        "summaries": {
            "wt": {
                "min": float(wt.min()),
                "max": float(wt.max()),
                "mean": float(wt.mean()),
                "sd": float(wt.std()),
            },
            "deletion": {
                "min": float(deletion.min()),
                "max": float(deletion.max()),
                "mean": float(deletion.mean()),
                "sd": float(deletion.std()),
            },
            "difference": {
                "min": float((deletion - wt).min()),
                "max": float((deletion - wt).max()),
                "mean": float((deletion - wt).mean()),
                "sd": float((deletion - wt).std()),
            },
        },
        "common_inputs_sha256": sha256(args.common_inputs),
    }
    (args.out_dir / "CHIMAERA_RUN_AUDIT.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit["summaries"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
