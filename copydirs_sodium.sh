#!/bin/bash


remote_base="/Volumes/DRS-Sodium-MRI/DATA/NaSCAR/analysis_output/PIPELINE"
local_base="/Volumes/nemosine/SAN/NASCAR"

for subj in Subject{1..4}; do
  for site in site{1..3}; do
    #src_dir="${remote_base}/${subj}/${site}/pipeline/reference_sodium"
    #dst_dir="${local_base}/${subj}/${site}/pipeline/reference_sodium"
    src_dir="${remote_base}/${subj}/${site}/pipeline/other_sodium"
    dst_dir="${local_base}/${subj}/${site}/pipeline/other_sodium"


    if [ -d "$src_dir" ]; then
      echo "🔍 Checking $src_dir"
      rsync -avz \
        --include='*TSC_3_bottles_mode.nii' \
        --exclude='*' \
        "$src_dir/" "$dst_dir/"

    # if [ -d "$src_dir" ]; then
    #   echo "🔍 Checking $src_dir"
    #   rsync -avz \
    #     --include='*TSC_3_bottles.nii' \
    #     --include='*TSC_B1.nii' \
    #     --exclude='*' \
    #     "$src_dir/" "$dst_dir/"
    else
      echo "⚠️ Missing remote folder: $src_dir"
    fi
  done
done
