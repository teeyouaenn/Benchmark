#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:?usage: download_resources.sh OUTPUT_ROOT}"
PYTHON="${PYTHON:-python}"
mkdir -p "$ROOT"/{akita_v2,deepc,epcot,chimaera}

fetch() {
  local url="$1"
  local out="$2"
  if [[ ! -s "$out" ]]; then
    curl --fail --location --retry 5 --retry-delay 5 "$url" --output "$out"
  fi
}

# AkitaV2: fold 2 is the held-out fold containing the BHLHE40 locus.
fetch "https://storage.googleapis.com/basenji_hic/3-2021/models/f2c0/train/model0_best.h5" \
  "$ROOT/akita_v2/model0_best.h5"
fetch "https://storage.googleapis.com/basenji_hic/3-2021/models/f2c0/train/params.json" \
  "$ROOT/akita_v2/params.json"
fetch "https://storage.googleapis.com/basenji_hic/3-2021/data/hg38/targets.txt" \
  "$ROOT/akita_v2/targets.txt"
fetch "https://storage.googleapis.com/basenji_hic/3-2021/data/hg38/sequences.bed" \
  "$ROOT/akita_v2/sequences.bed"

# DeepC: highest-resolution official K562 checkpoint.
fetch "https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/models/model_deepCregr_5kb_K562.tar.gz" \
  "$ROOT/deepc/model_deepCregr_5kb_K562.tar.gz"
if [[ ! -d "$ROOT/deepc/model_deepCregr_5kb_K562" ]]; then
  tar -xzf "$ROOT/deepc/model_deepCregr_5kb_K562.tar.gz" -C "$ROOT/deepc"
fi

# EPCOT: 1-kb HFF Micro-C head plus the released pretraining checkpoint.
if [[ ! -s "$ROOT/epcot/HFF_Micro-C_transformer.pt" ]]; then
  "$PYTHON" -m gdown 1PUQyBdqadq2AI9IZPpeAYc5zwjWeoph7 \
    -O "$ROOT/epcot/HFF_Micro-C_transformer.pt"
fi
if [[ ! -s "$ROOT/epcot/pretrain_dnase.pt" ]]; then
  "$PYTHON" -m gdown 1_YfpNSv-2ABQV2qSyBxem-y7aJFyRNzz \
    -O "$ROOT/epcot/pretrain_dnase.pt"
fi

# Chimaera: official human model/data release linked from the official Colab.
fetch "https://osf.io/cwpke/download" "$ROOT/chimaera/human_release.tar.gz"

find "$ROOT" -type f -print0 | sort -z | xargs -0 sha256sum > "$ROOT/SHA256SUMS.txt"
