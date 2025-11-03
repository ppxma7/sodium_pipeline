#!/bin/bash
# Wrapper for run_sodium_pipeline.py
# Usage: ./run_sodium_pipeline.sh SUBJECT_ID

if [ $# -ne 1 ]; then
    echo "Usage: $0 SUBJECT_ID"
    echo "Example: $0 16905_004"
    exit 1
fi

SUBJECT=$1

# Root directory where subject folders live
ROOTDIR="/Volumes/nemosine/SAN/SASHB/inputs"

# Build paths
ARG1="${ROOTDIR}/${SUBJECT}/MPRAGE"
ARG2="${ROOTDIR}/${SUBJECT}/MPRAGE_sodium"
ARG3="${ROOTDIR}/${SUBJECT}/sodium"

# Check that all three paths exist
for d in "$ARG1" "$ARG2" "$ARG3"; do
    if [ ! -d "$d" ]; then
        echo "Error: directory $d does not exist."
        exit 1
    fi
done

echo "Running sodium pipeline for $SUBJECT"
python3 /Users/ppzma/Documents/MATLAB/mycode/sodium_pipeline/run_sodium_pipeline.py "$ARG1" "$ARG2" "$ARG3"
