import os
import subprocess
import argparse
import sys
import shutil
import glob
import nibabel as nib
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET


# ---------- USER CONFIG ----------
FSLDIR = "/usr/local/fsl"
MNI_TEMPLATE = f"{FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz"
MNI_BRAIN_MASK = f"{FSLDIR}/data/standard/MNI152_T1_1mm_brain_mask.nii.gz"
MY_CONFIG_DIR = "/Users/ppzma/data"  # contains bb_fnirt.cnf

parser = argparse.ArgumentParser(description="Run sodium MRI pipeline")

parser.add_argument("ARG1", help="Path to outputs")
parser.add_argument("ARG2", help="Path to reference sodium")
parser.add_argument("ARG3", help="Path to proton images")


args = parser.parse_args()
ARG1 = args.ARG1
ARG2 = args.ARG2
ARG3 = args.ARG3

subject = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(ARG1))))
print(subject)


# --- Auto-detect input files ---
# 1. Main MPRAGE (should start with WIP_MPRAGE_ or similar)
mprage_matches = glob.glob(os.path.join(ARG1, "*MPRAGE_optibrain.nii*"))
if len(mprage_matches) == 0:
    raise FileNotFoundError(f"No MPRAGE file found in {ARG1}")
elif len(mprage_matches) > 1:
    print(f"⚠️ Multiple MPRAGE files found, using first: {mprage_matches[0]}")
mprage_file = mprage_matches[0]

sodium_matches = glob.glob(os.path.join(ARG1, "*align12dof.nii.gz"))
if len(sodium_matches) == 0:
    raise FileNotFoundError(f"No sodium align12dof files found in {ARG1}")
else:
	sodium_files = sorted(sodium_matches)
	print("Found sodium files:", sodium_matches)
	print(f"Found {len(sodium_matches)} sodium files")


# mprage paths
affine_mprage_to_mni = os.path.join(ARG1, f"{subject}_mprage2mni.mat")
mprage_to_mni = os.path.join(ARG1, f"{subject}_mprage2mni_linear.nii.gz")


def strip_ext(fname):
    if fname.endswith(".nii.gz"):
        return fname[:-7]  # strip 7 chars for '.nii.gz'
    else:
        return os.path.splitext(fname)[0]


def run(cmd, check=True):
    print("🔧 Running:", " ".join(cmd))
    subprocess.run(cmd, check=check)


# Pick the first sodium file as the "moving" image
#ref_sodium = sodium_files[0]
#matches = glob.glob(os.path.join(ARG1, f"{ARG2}.nii*"))
matches = []
for pattern in [f"{ARG2}.nii*", f"{ARG2}_2375.nii*"]:
    found = glob.glob(os.path.join(ARG1, pattern))
    matches.extend(found)

if not matches:
    raise FileNotFoundError(f"Reference sodium {ARG2}.nii* or {ARG2}_2375.nii* not found in {ARG1}")

# Prefer resampled if both exist
if any("_2375" in m for m in matches):
    ref_sodium = [m for m in matches if "_2375" in m][0]
else:
    ref_sodium = matches[0]

print(f"✅ Using reference sodium file: {ref_sodium}")

#ref_sodium = matches[0]


# --- Find matching TSC file ---
ref_base = strip_ext(ref_sodium)
# Handle both “base_TSC” and “base_2375_TSC” forms
tsc_patterns = [
    f"{ref_base}_TSC.nii*",
    f"{ref_base.replace('_2375', '')}_TSC_2375.nii*"
]
tsc_matches = []
for pat in tsc_patterns:
    tsc_matches.extend(glob.glob(pat))

if not tsc_matches:
    print(f"ℹ️ No TSC file found for {ref_sodium}, continuing without TSC")
    ref_sodium_tsc = None
else:
    # Prefer a matching _2375 version if it exists
    if any("_2375" in t for t in tsc_matches):
        ref_sodium_tsc = [t for t in tsc_matches if "_2375" in t][0]
    else:
        ref_sodium_tsc = tsc_matches[0]
    print(f"✅ Found reference TSC file: {ref_sodium_tsc}")

    

# Look for reference sodium TSC
# ref_base = strip_ext(ref_sodium)
# tsc_matches = glob.glob(f"{ref_base}_TSC.nii*")
# if not tsc_matches:
#     print(f"ℹ️ No TSC file found for {ref_sodium}, continuing without TSC")
#     ref_sodium_tsc = None
# else:
#     ref_sodium_tsc = tsc_matches[0]
#     print(f"✅ Found reference TSC file: {ref_sodium_tsc}")



base = strip_ext(ref_sodium)
sodium_2_sodiumMPRAGE = f"{base}_toMPRAGE.nii.gz"
sodium_2_sodiumMPRAGE_mat = f"{base}_toMPRAGE.mat"

#####
# First lets move sodium files to MPRAGE



if not os.path.exists(sodium_2_sodiumMPRAGE):
    run([
        f"{FSLDIR}/bin/flirt",
        "-in", ref_sodium,
        "-ref", mprage_file,
        "-omat", sodium_2_sodiumMPRAGE_mat,
        "-out", sodium_2_sodiumMPRAGE,
        "-bins", "256",
        "-dof", "6",
        "-schedule", f"{FSLDIR}/etc/flirtsch/sch3Dtrans_3dof",
        "-cost", "normmi",
        "-searchrx", "0", "0",
        "-searchry", "0", "0",
        "-searchrz", "0", "0",
        "-interp", "trilinear"
    ])
    print(f"✅ Ran FLIRT for {ref_sodium} → {sodium_2_sodiumMPRAGE}")
else:
    print("⏭️ Sodium→MPRAGE already exists, skipping FLIRT.")

#FOR SITE 3, using 12DOF for some reason
# if not os.path.exists(sodium_2_sodiumMPRAGE):
#     run([
#         f"{FSLDIR}/bin/flirt",
#         "-in", ref_sodium,
#         "-ref", mprage_file,
#         "-omat", sodium_2_sodiumMPRAGE_mat,
#         "-out", sodium_2_sodiumMPRAGE,
#         "-bins", "256",
#         "-dof", "12",
#         "-cost", "normmi",
#         "-searchrx", "0", "0",
#         "-searchry", "0", "0",
#         "-searchrz", "0", "0",
#         "-interp", "trilinear"
#     ])
#     print(f"✅ Ran FLIRT for {ref_sodium} → {sodium_2_sodiumMPRAGE}")
# else:
#     print("⏭️ Sodium→MPRAGE already exists, skipping FLIRT.")


#sys.exit(0)

# Apply the same transform to all the *other* sodium files
for s in sodium_files:
    out_file = f"{strip_ext(s)}_toMPRAGE.nii.gz"
    if os.path.exists(out_file):
        print(f"⏭️ Skipping {os.path.basename(s)} (already aligned).")
        continue
    run([
        f"{FSLDIR}/bin/flirt",
        "-in", s,
        "-ref", mprage_file,
        "-applyxfm",
        "-init", sodium_2_sodiumMPRAGE_mat,
        "-out", out_file
    ])
    print(f"✅ Applied transform to {os.path.basename(s)} → {out_file}")

if ref_sodium_tsc:
    base_tsc = strip_ext(ref_sodium_tsc)
    ref_sodium_tsc_toMPRAGE = f"{base_tsc}_toMPRAGE.nii.gz"

    if not os.path.exists(ref_sodium_tsc_toMPRAGE):
        run([
            f"{FSLDIR}/bin/flirt",
            "-in", ref_sodium_tsc,
            "-ref", mprage_file,
            "-applyxfm",
            "-init", sodium_2_sodiumMPRAGE_mat,  # use same transform as sodium
            "-out", ref_sodium_tsc_toMPRAGE
        ])
        print(f"✅ Registered reference TSC to MPRAGE: {ref_sodium_tsc_toMPRAGE}")


# MPRAGE to MNI space
if not os.path.exists(mprage_to_mni):
    print("Running FLIRT MPRAGE → MNI")
    run([
        f"{FSLDIR}/bin/flirt",
        "-in", mprage_file,
        "-ref", MNI_TEMPLATE,
        "-omat", affine_mprage_to_mni,
        "-out", mprage_to_mni,
        "-bins", "256",
        "-dof", "12",
        "-cost", "corratio",
        "-searchrx", "-90", "90",
        "-searchry", "-90", "90",
        "-searchrz", "-90", "90",
        "-interp", "trilinear"
    ])
else:
    print("⏭️ FLIRT MPRAGE→MNI exists, skipping.")

# Collect all sodium images already aligned to MPRAGE:
# (reference sodium + the others)
all_sodium_in_mprage = [sodium_2_sodiumMPRAGE] + [
    f"{strip_ext(s)}_toMPRAGE.nii.gz" for s in sodium_files
]

# Add reference TSC if it was aligned
if ref_sodium_tsc:
    base_tsc = strip_ext(ref_sodium_tsc)
    ref_sodium_tsc_toMPRAGE = f"{base_tsc}_toMPRAGE.nii.gz"
    all_sodium_in_mprage.append(ref_sodium_tsc_toMPRAGE)


# Apply MPRAGE→MNI to each
for sodium_in_mprage_space in all_sodium_in_mprage:
    if not os.path.exists(sodium_in_mprage_space):
        print(f"⚠️ Missing MPRAGE-aligned sodium: {sodium_in_mprage_space}, skipping.")
        continue

    sodium_file_mni = f"{strip_ext(sodium_in_mprage_space)}_MNI.nii.gz"

    if os.path.exists(sodium_file_mni):
        print(f"⏭️ Skipping {os.path.basename(sodium_in_mprage_space)} (already in MNI space).")
        continue

    run([
        f"{FSLDIR}/bin/flirt",
        "-in", sodium_in_mprage_space,
        "-ref", MNI_TEMPLATE,
        "-applyxfm",
        "-init", affine_mprage_to_mni,
        "-out", sodium_file_mni
    ])
    print(f"✅ Sodium moved to MNI space: {sodium_file_mni}")


#### MOVING Place outputs_mni one level up from ARG1
parent_dir = os.path.dirname(os.path.dirname(ARG1))
outputs_mni = os.path.join(parent_dir, "outputs_mni")
os.makedirs(outputs_mni, exist_ok=True)

mni_files = glob.glob(os.path.join(ARG1, "*MNI*.nii*"))

for f in mni_files:
    if not os.path.exists(f):
        print(f"⚠️ Missing source file: {f}")
        continue

    dest = os.path.join(outputs_mni, os.path.basename(f))
    if os.path.exists(dest):
        print(f"⏭️ Already exists in outputs_mni: {dest}")
        continue

    shutil.copy(f, dest)
    print(f"📦 Copied {os.path.basename(f)} → {outputs_mni}/")


mprage_mni_file = os.path.join(ARG1, f"{subject}_mprage2mni_linear.nii.gz")
dest_file = os.path.join(outputs_mni, os.path.basename(mprage_mni_file))

if not os.path.exists(mprage_mni_file):
    print(f"⚠️ Missing file: {mprage_mni_file}")
elif os.path.exists(dest_file):
    print(f"⏭️ Already copied: {dest_file}")
else:
    shutil.copy(mprage_mni_file, dest_file)
    print(f"📦 Copied {os.path.basename(mprage_mni_file)} → {outputs_mni}/")

##### move atlas to sodium

atlas = f"{FSLDIR}/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz"
    

# Step 1: Invert MPRAGE→MNI (get MNI→MPRAGE)
mni2mprage_mat = os.path.join(ARG1, f"{subject}_mni2mprage_linear.mat")
if not os.path.exists(mni2mprage_mat):
    run([
        f"{FSLDIR}/bin/convert_xfm",
        "-omat", mni2mprage_mat,
        "-inverse", affine_mprage_to_mni
    ])
    print(f"✅ Created MNI→MPRAGE matrix: {mni2mprage_mat}")
else:
    print("⏭️ MNI→MPRAGE matrix already exists.")

# Step 2: Invert SODIUM→SODIUMMPRAGE (get SODIUMMPRAGE→SODIUM)
mprage2sodium_mat = os.path.join(ARG1, f"{subject}_mprage2sodium.mat")
if not os.path.exists(mprage2sodium_mat):
    run([
        f"{FSLDIR}/bin/convert_xfm",
        "-omat", mprage2sodium_mat,
        "-inverse", sodium_2_sodiumMPRAGE_mat
    ])
    print(f"✅ Created MPRAGE→Sodium matrix: {mprage2sodium_mat}")

#step 3
mni2sodium_mat = os.path.join(ARG1, f"{subject}_mni2sodium.mat")
if not os.path.exists(mni2sodium_mat):
    run([
        f"{FSLDIR}/bin/convert_xfm",
        "-omat", mni2sodium_mat,
        "-concat", mprage2sodium_mat,
        mni2mprage_mat
    ])
    print(f"✅ Created MNI→Sodium matrix: {mni2sodium_mat}")


# Step : Apply transform to atlas (nearest neighbour for labels)
atlas_in_sodium = os.path.join(ARG1, f"{subject}_atlas_in_sodium.nii.gz")


if not os.path.exists(atlas_in_sodium):
    run([
        f"{FSLDIR}/bin/flirt",
        "-in", atlas,
        "-ref", ref_sodium,   # use the sodium reference image as target
        "-applyxfm",
        "-init", mni2sodium_mat,
        "-interp", "nearestneighbour",
        "-out", atlas_in_sodium
    ])
    print(f"✅ Atlas moved to sodium space: {atlas_in_sodium}")

# Also want to move the FAST outputs to sodium space, so that's just applying 
# the MPRAGE->Sodium matrix to the FAST outputs
# Proton images are in ARG3

# Collect all FAST outputs
fast_matches = []
fast_matches += glob.glob(os.path.join(ARG3, "*MPRAGE_optibrain_pve*.nii*"))
fast_matches += glob.glob(os.path.join(ARG3, "*MPRAGE_optibrain_seg.nii*"))
fast_matches += glob.glob(os.path.join(ARG3, "*MPRAGE_optibrain_mixeltype.nii*"))

# Loop over each file and apply the transform
for fast_file in fast_matches:
    base = os.path.basename(fast_file)
    out_file = os.path.join(ARG3, base.replace("MPRAGE_optibrain", "fast_in_sodium"))

    if not os.path.exists(out_file):
        run([
            f"{FSLDIR}/bin/flirt",
            "-in", fast_file,
            "-ref", ref_sodium,
            "-applyxfm",
            "-init", mprage2sodium_mat,
            "-interp", "nearestneighbour",
            "-out", out_file
        ])
        print(f"✅ Transformed {base} → {out_file}")
    else:
        print(f"Skipping {out_file} (already exists)")

# --- Binarise the pve masks and apply them to TSC files ---
print("\n--- Processing FAST masks and applying to TSC files ---")

# 1. Get the PVE mask files we just generated in sodium space (but exclude _bin files)
pve_masks = sorted(
    f for f in glob.glob(os.path.join(ARG3, "*fast_in_sodium_pve_*.nii.gz"))
    if "_bin" not in f
)
if not pve_masks:
    raise FileNotFoundError(f"No PVE masks found in {ARG3}")


# Optional: define a binarisation threshold (typical ~0.5)
PVE_THRESHOLD = 0.5

# 2. Binarise each PVE mask
binarised_masks = []
for pve_file in pve_masks:
    bin_file = pve_file.replace(".nii.gz", "_bin.nii.gz")
    if not os.path.exists(bin_file):
        run([
            f"{FSLDIR}/bin/fslmaths",
            pve_file,
            "-thr", str(PVE_THRESHOLD),
            "-bin",
            bin_file
        ])
        print(f"✅ Binarised {pve_file} → {bin_file}")
    else:
        print(f"Skipping {bin_file} (already exists)")
    binarised_masks.append(bin_file)


# --- 2a. Erode the PVE2 masks (WM) ---
print("\n--- Eroding PVE2 (WM) masks ---")
pve2_files = [f for f in binarised_masks if "_pve_2_bin.nii.gz" in f]
for pve2_file in pve2_files:
    pve2_ero_file = pve2_file.replace(".nii.gz", "_ero.nii.gz")
    if not os.path.exists(pve2_ero_file):
        run([
            f"{FSLDIR}/bin/fslmaths",
            pve2_file,
            "-kernel", "sphere", "3.5",
            "-ero",
            pve2_ero_file
        ])
        print(f"✅ Eroded {pve2_file} → {pve2_ero_file}")
    else:
        print(f"Skipping {pve2_ero_file} (already exists)")

    # Replace the non-eroded file in binarised_masks
    idx = binarised_masks.index(pve2_file)
    binarised_masks[idx] = pve2_ero_file


# 3. Find TSC files to apply masks to
#all_tsc_files = glob.glob(os.path.join(ARG1, "*TSC*.nii.gz"))
#tsc_files = [f for f in all_tsc_files if "MPRAGE" not in f and "MNI" not in f]
all_tsc_files = [
    f for f in glob.glob(os.path.join(ARG1, "*TSC*.nii*"))
    if "_masked_" not in f   # ✅ don't mask masked files again
]

tsc_files = [
    f for f in all_tsc_files
    if "MPRAGE" not in f
    and "MNI" not in f
    and "bet" not in f
    and (
        ("align" in f or "resampled" in f)
        or os.path.basename(f).endswith("_TSC.nii")
        or os.path.basename(f).endswith("_TSC.nii.gz")
    )
]

# tsc_files = [
#     f for f in all_tsc_files
#     if "MPRAGE" not in f
#     and "MNI" not in f
#     and "bet" not in f
#     and ("align" in f or "resampled" in f)
# ]

if not tsc_files:
    print("⚠️ No matching TSC files found.")
else:
    print(f"Found {len(tsc_files)} TSC files to mask.")



# 4. Apply each binary mask to each TSC file
for tsc_file in tsc_files:
    tsc_base = strip_ext(os.path.basename(tsc_file))
    tsc_img = nib.load(tsc_file)
    tsc_shape = tsc_img.shape

    for mask in binarised_masks:
        mask_base = os.path.basename(mask).replace(".nii.gz", "")
        # Strip any trailing _bin to avoid _bin_bin_bin outputs

        mask_base = mask_base.replace("_bin", "").replace("_ero", "")
        out_file = os.path.join(ARG1, f"{tsc_base}_masked_{mask_base}_bin.nii.gz")

        mask_img = nib.load(mask)
        mask_shape = mask_img.shape

        # --- Handle shape mismatch ---
        if tsc_shape != mask_shape:
            print(f"⚠️ Shape mismatch for {tsc_base}: TSC {tsc_shape} vs Mask {mask_shape} — skipping")
            continue

        if not os.path.exists(out_file):
            run([
                f"{FSLDIR}/bin/fslmaths",
                tsc_file,
                "-mas", mask,
                out_file
            ])
            print(f"✅ Masked {tsc_file} with {mask_base} → {out_file}")
        else:
            print(f"Skipping {out_file} (already exists)")


#### MOVING Place outputs_mni one level up from ARG1
#parent_dir = os.path.dirname(os.path.dirname(ARG1))
outputs_pve_native = os.path.join(parent_dir, "outputs_pve_native")
os.makedirs(outputs_pve_native, exist_ok=True)

pve_files = glob.glob(os.path.join(ARG1, "*_pve_*.nii*"))

for f in pve_files:
    if not os.path.exists(f):
        print(f"⚠️ Missing source file: {f}")
        continue

    dest = os.path.join(outputs_pve_native, os.path.basename(f))
    # if os.path.exists(dest):
    #     print(f"⏭️ Already exists in outputs_pve_native: {dest}")
    #     continue

    shutil.move(f, dest)
    print(f"📦 Moved {os.path.basename(f)} → {outputs_pve_native}/")


# --- Move FAST outputs to MNI space as well
print("\n--- Transforming FAST outputs to MNI space ---")

for fast_file in fast_matches:
    base = os.path.basename(fast_file)
    out_file = os.path.join(ARG3, base.replace("MPRAGE_optibrain", "fast_in_MNI"))

    if not os.path.exists(out_file):
        run([
            f"{FSLDIR}/bin/flirt",
            "-in", fast_file,
            "-ref", MNI_TEMPLATE,
            "-applyxfm",
            "-init", affine_mprage_to_mni,
            "-interp", "nearestneighbour",
            "-out", out_file
        ])
        print(f"✅ Transformed {base} → {out_file}")
    else:
        print(f"⏭️ Skipping {out_file} (already exists)")

# --- Binarise the pve masks in MNI space ---
print("\n--- Processing FAST masks ---")

# 1. 
pve_masks = sorted(
    f for f in glob.glob(os.path.join(ARG3, "*fast_in_MNI_pve_*.nii.gz"))
    if "_bin" not in f
)
if not pve_masks:
    raise FileNotFoundError(f"No PVE masks found in {ARG3}")

# 2. Binarise each PVE mask
binarised_masks = []
for pve_file in pve_masks:
    bin_file = pve_file.replace(".nii.gz", "_bin.nii.gz")
    if not os.path.exists(bin_file):
        run([
            f"{FSLDIR}/bin/fslmaths",
            pve_file,
            "-thr", str(PVE_THRESHOLD),
            "-bin",
            bin_file
        ])
        print(f"✅ Binarised {pve_file} → {bin_file}")
    else:
        print(f"Skipping {bin_file} (already exists)")
    binarised_masks.append(bin_file)


# --- 2a. Erode the PVE2 masks (WM) ---
print("\n--- Eroding PVE2 (WM) masks ---")
pve2_files = [f for f in binarised_masks if "_pve_2_bin.nii.gz" in f]
for pve2_file in pve2_files:
    pve2_ero_file = pve2_file.replace(".nii.gz", "_ero.nii.gz")
    if not os.path.exists(pve2_ero_file):
        run([
            f"{FSLDIR}/bin/fslmaths",
            pve2_file,
            "-kernel", "sphere", "2",
            "-ero",
            pve2_ero_file
        ])
        print(f"✅ Eroded {pve2_file} → {pve2_ero_file}")
    else:
        print(f"Skipping {pve2_ero_file} (already exists)")

    # Replace the non-eroded file in binarised_masks
    idx = binarised_masks.index(pve2_file)
    binarised_masks[idx] = pve2_ero_file



#### MOVING MNI FAST OUTPUTS
outputs_pve_mni = os.path.join(parent_dir, "outputs_pve_mni")
os.makedirs(outputs_pve_mni, exist_ok=True)

pve_files = glob.glob(os.path.join(ARG3, "*MNI_pve_*.nii*"))

for f in pve_files:
    if not os.path.exists(f):
        print(f"⚠️ Missing source file: {f}")
        continue
    dest = os.path.join(outputs_pve_mni, os.path.basename(f))
    shutil.copy(f, dest)
    print(f"📦 Copied {os.path.basename(f)} → {outputs_pve_mni}/")

