#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "usage: $0 REPOSITORY CHECKPOINT_PREFIX FASTA OUT_DIR AUDIT_SCRIPT" >&2
  exit 2
fi

repository=$1
checkpoint=$2
fasta=$3
out_dir=$4
audit_script=$5
python_bin="${DEEPC_PYTHON:-python}"
source_dir="$repository/tensorflow2.0_compatibility_version"

mkdir -p "$out_dir"
printf 'chr3\t4976067\t4976790\treference\n' > "$out_dir/wt.variant.tsv"
printf 'chr3\t4976067\t4976790\t.\n' > "$out_dir/deletion.variant.tsv"

site=$($python_bin -c 'import site; print(site.getsitepackages()[0])')
cuda_libs=$(find "$site/nvidia" -type d -name lib 2>/dev/null | paste -sd:)
if [[ -n "$cuda_libs" ]]; then
  export LD_LIBRARY_PATH="${cuda_libs}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi
export TF_FORCE_GPU_ALLOW_GROWTH=true

cd "$source_dir"
for condition in wt deletion; do
  CUDA_VISIBLE_DEVICES=0 "$python_bin" run_deploy_shape_deepCregr.py \
    --dlmodel deepCregr \
    --batch_size 1 \
    --out_dir "$out_dir/$condition" \
    --name_tag "$condition" \
    --input "$out_dir/$condition.variant.tsv" \
    --model "$checkpoint" \
    --genome "$fasta" \
    --bp_context 1005000 \
    --add_window 200000 \
    --num_classes 201 \
    --bin_size 5000 \
    --bin_steps full \
    --run_on gpu \
    --gpu 0
done

"$python_bin" "$audit_script" --out-dir "$out_dir"
