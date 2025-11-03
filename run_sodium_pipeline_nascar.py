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
OPTIBET_PATH = "/Users/ppzma/Documents/MATLAB/optibet.sh"
FASTPATH = f"{FSLDIR}/share/fsl/bin/fast"

# ---------- ARGPARSE ----------
parser = argparse.ArgumentParser(description="Run sodium MRI pipeline")
parser.add_argument("ARG1", help="Path to MPRAGE directory")
parser.add_argument("ARG2", help="Path to reference sodium")
parser.add_argument("ARG3", help="Path to other sodium")
parser.add_argument("ARG4", help="Basename of reference sodium (e.g. floret, radial, seiffert)")

args = parser.parse_args()

ARG1 = args.ARG1
ARG2 = args.ARG2
ARG3 = args.ARG3
ARG4 = args.ARG4

def get_extension(fname):
    if fname.endswith(".nii.gz"):
        return ".nii.gz"
    return os.path.splitext(fname)[1]


subject = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(ARG1))))
print(f"Running: {subject}")

# --- Auto-detect input files ---
# 1. Main MPRAGE (should start with WIP_MPRAGE_ or similar)
mprage_matches = glob.glob(os.path.join(ARG1, "*3DT1*.nii*"))
if len(mprage_matches) == 0:
    raise FileNotFoundError(f"No MPRAGE file found in {ARG1}")
elif len(mprage_matches) > 1:
    print(f"Multiple MPRAGE files found, using first: {mprage_matches[0]}")
mprage_file = mprage_matches[0]

print(f"This is the T1w: {mprage_file}")

# 2. Reference sodium
sodium_ref_matches = glob.glob(os.path.join(ARG2, f"{ARG4}.nii*"))
if len(sodium_ref_matches) == 0:
    raise FileNotFoundError(f"No sodium ref {ARG4}.nii found in {ARG2}")
elif len(sodium_ref_matches) > 1:
    print(f"Multiple sodium ref files found, using first: {sodium_ref_matches[0]}")
sodium_ref_file = sodium_ref_matches[0]


# 3. Reference sodium TSC
ref_base = os.path.splitext(os.path.basename(sodium_ref_file))[0]
sodium_ref_tsc_matches = glob.glob(os.path.join(ARG2, f"{ref_base}_TSC.nii*"))
if len(sodium_ref_tsc_matches) == 0:
    print(f"No sodium ref TSC file found in {ARG2} for {ref_base}, continuing without TSC")
    sodium_ref_tsc_file = None
elif len(sodium_ref_tsc_matches) > 1:
    print(f"Multiple sodium ref TSC files found, using first: {sodium_ref_tsc_matches[0]}")
    sodium_ref_tsc_file = sodium_ref_tsc_matches[0]
else:
    sodium_ref_tsc_file = sodium_ref_tsc_matches[0]


#print(sodium_ref_file)
print(f"This is the sodium reference: {sodium_ref_file}")

#print(sodium_ref_tsc_file)

# 4. Grab the other sodium files
sodium_types = ["floret", "radial", "seiffert"]
# Remove the chosen reference to get the "other two"
others = [s for s in sodium_types if s.lower() != ARG4.lower()]
other_sodium_files = {}


for oth in others:
    if oth == "seiffert":
        # Seiffert 
        matches = glob.glob(os.path.join(ARG3, "*23NaSeiffert*.nii*"))
    elif oth == "radial":
        # Radial 
        #matches = glob.glob(os.path.join(ARG3, "radialfiltstrong.nii*"))
        matches = []
        for pattern in ["radialfiltstrong.nii*", "radialdefault.nii*"]:
            matches.extend(glob.glob(os.path.join(ARG3, pattern)))
    elif oth == "floret":
        # Floret 
        matches = glob.glob(os.path.join(ARG3, "floret.nii*"))
    else:
        matches = []

    if len(matches) == 0:
        print(f"Could not find expected sodium file for {oth} in {ARG3}")
        other_sodium_files[oth] = None
        continue

   
    if len(matches) > 1:
        print(f"Multiple matches for {oth}, using first: {matches[0]}")

    src = matches[0]

    # Copy radial and seiffert to standardized names
    if oth in ["seiffert", "radial", "floret"]:
        ext = get_extension(src)  # e.g. ".nii.gz"
        dst = os.path.join(ARG3, f"{oth}{ext}")  # keep proper extension
        #dst = os.path.join(ARG3, f"{oth}.nii*")
        if os.path.exists(dst):
            print(f"⏭️ Skipping copy for {oth}, {dst} already exists")
        else:
            shutil.copy(src, dst)
            print(f"✅ Copied {oth}: {src} → {dst}")
        other_sodium_files[oth] = dst
    else:
        # Keep floret as-is
        print(f"✅ Found {oth}: {src}")
        #other_sodium_files[oth] = src

# Filter out missing files (None values)
valid_other_sodiums = {k: v for k, v in other_sodium_files.items() if v is not None}
print(valid_other_sodiums)



# 5. Reference sodium TSC
ref_base = os.path.splitext(os.path.basename(sodium_ref_file))[0]
sodium_ref_tsc_matches = glob.glob(os.path.join(ARG2, f"{ref_base}_TSC.nii*"))
if len(sodium_ref_tsc_matches) == 0:
    print(f"No sodium ref TSC file found in {ARG2} for {ref_base}, continuing without TSC")
    sodium_ref_tsc_file = None
elif len(sodium_ref_tsc_matches) > 1:
    print(f"Multiple sodium ref TSC files found, using first: {sodium_ref_tsc_matches[0]}")
    sodium_ref_tsc_file = sodium_ref_tsc_matches[0]
else:
    sodium_ref_tsc_file = sodium_ref_tsc_matches[0]

# --- Find TSC files for the two "other sodium" images ---
other_tsc_files = {}
for name, f in valid_other_sodiums.items():
    base = f[:-7] if f.endswith(".nii.gz") else os.path.splitext(f)[0]
    tsc_matches = glob.glob(f"{base}_TSC.nii*")

    if len(tsc_matches) == 0:
        print(f"No TSC file found for {name} ({f}), continuing without TSC")
        other_tsc_files[name] = None
    elif len(tsc_matches) > 1:
        print(f"Multiple TSC files for {name}, using first: {tsc_matches[0]}")
        other_tsc_files[name] = tsc_matches[0]
    else:
        other_tsc_files[name] = tsc_matches[0]

#print(other_tsc_files)

#sys.exit(0)

######################
# Okay here we have all our files. the sodium images (file1 and file2)
# The reference sodium image
# the 3dt1

######################


def run(cmd, check=True):
    print("🔧 Running:", " ".join(cmd))
    subprocess.run(cmd, check=check)

def strip_ext(fname):
    if fname.endswith(".nii.gz"):
        return fname[:-7]  # strip 7 chars for '.nii.gz'
    else:
        return os.path.splitext(fname)[0]



# # Step 1) Resample Seiffert if present
# if other_sodium_files.get("seiffert"):
#     seiffert_file = other_sodium_files["seiffert"]
#     resampled_seiffert = f"{os.path.splitext(seiffert_file)[0]}_2375.nii.gz"
#     print(resampled_seiffert)    
#     if os.path.exists(resampled_seiffert):
#         print(f" Skipping Seiffert resample, already exists: {resampled_seiffert}")
#     else:
#         run(["flirt", "-in", seiffert_file, "-ref", seiffert_file,
#              "-out", resampled_seiffert, "-applyisoxfm", "2.375"])

#     if other_tsc_files.get("seiffert"):
#         seiffert_tsc = other_tsc_files["seiffert"]
#         resampled_seiffert_tsc = f"{os.path.splitext(seiffert_tsc)[0]}_2375.nii.gz"

#         if os.path.exists(resampled_seiffert_tsc):
#             print(f" Skipping Seiffert TSC resample, already exists: {resampled_seiffert_tsc}")
#         else:
#             run(["flirt", "-in", seiffert_tsc, "-ref", seiffert_tsc,
#                  "-out", resampled_seiffert_tsc, "-applyisoxfm", "2.375"])
# else:
#     print(" No Seiffert sodium file found, skipping resampling.")


# Step 1) Resample Seiffert if present (reference or other)
resampled_seiffert = None
resampled_seiffert_tsc = None

# Case 1: Seiffert is the reference sodium
if ARG4.lower() == "seiffert":
    print("Seiffert is the ref")
    seiffert_file = sodium_ref_file
    resampled_seiffert = f"{strip_ext(seiffert_file)}_2375.nii.gz"

    if os.path.exists(resampled_seiffert):
        print(f"⏭️ Skipping Seiffert resample (reference), already exists: {resampled_seiffert}")
    else:
        run(["flirt", "-in", seiffert_file, "-ref", seiffert_file,
             "-out", resampled_seiffert, "-applyisoxfm", "2.375"])

    if sodium_ref_tsc_file:
        seiffert_tsc = sodium_ref_tsc_file
        resampled_seiffert_tsc = f"{strip_ext(seiffert_tsc)}_2375.nii.gz"
        if os.path.exists(resampled_seiffert_tsc):
            print(f"⏭️ Skipping Seiffert TSC resample (reference), already exists: {resampled_seiffert_tsc}")
        else:
            run(["flirt", "-in", seiffert_tsc, "-ref", seiffert_tsc,
                 "-out", resampled_seiffert_tsc, "-applyisoxfm", "2.375"])
    # now we need to repoint our ref to the seiffert resampled
    sodium_ref_file = resampled_seiffert

# Case 2: Seiffert is one of the "other sodiums"
elif other_sodium_files.get("seiffert"):
    seiffert_file = other_sodium_files["seiffert"]
    resampled_seiffert = f"{strip_ext(seiffert_file)}_2375.nii.gz"

    if os.path.exists(resampled_seiffert):
        print(f"⏭️ Skipping Seiffert resample, already exists: {resampled_seiffert}")
    else:
        run(["flirt", "-in", seiffert_file, "-ref", seiffert_file,
             "-out", resampled_seiffert, "-applyisoxfm", "2.375"])

    if other_tsc_files.get("seiffert"):
        seiffert_tsc = other_tsc_files["seiffert"]
        resampled_seiffert_tsc = f"{strip_ext(seiffert_tsc)}_2375.nii.gz"
        if os.path.exists(resampled_seiffert_tsc):
            print(f"⏭️ Skipping Seiffert TSC resample, already exists: {resampled_seiffert_tsc}")
        else:
            run(["flirt", "-in", seiffert_tsc, "-ref", seiffert_tsc,
                 "-out", resampled_seiffert_tsc, "-applyisoxfm", "2.375"])

else:
    print("⏭️ No Seiffert sodium file found (neither reference nor other), skipping resampling.")


# Do a check here to grab the relevant files
# Problem here though is what if seiffert is reference
radial_file = other_sodium_files.get("radial")
floret_file = other_sodium_files.get("floret")
seiffert_other = other_sodium_files.get("seiffert")
#print(seiffert_other)

#sys.exit(0)

# if radial_file:
#     print(f"Radial file: {radial_file}")
# else:
#     print("ℹNo radial file (probably reference or missing)")

# if floret_file:
#     print(f"Floret file: {floret_file}")
# else:
#     print("ℹNo floret file (probably reference or missing)")
print(f"Reminder, this is the sodium ref: {sodium_ref_file}")
#sys.exit(0)

radial_tsc_file = other_tsc_files.get("radial")
floret_tsc_file = other_tsc_files.get("floret")
seiffert_other_tsc = other_tsc_files.get("seiffert")

# resampled_seiffert_tsc

#2) BET sodium datasets which not reference
if radial_file:
    radial_bet = f"{strip_ext(radial_file)}_bet.nii.gz"
    if os.path.exists(radial_bet):
        print(f"⏭️ Skipping radial_bet, already exists: {radial_bet}")
    else:
        run(["bet", radial_file, radial_bet, "-f", "0.7"])

if floret_file:
    floret_bet = f"{strip_ext(floret_file)}_bet.nii.gz"
    if os.path.exists(floret_bet):
        print(f"⏭️ Skipping floret_bet, already exists: {floret_bet}")
    else:
        run(["bet", floret_file, floret_bet, "-f", "0.7"])

if seiffert_other:
    seiffert_bet = f"{strip_ext(resampled_seiffert)}_bet.nii.gz"
    if os.path.exists(seiffert_bet):
        print(f"⏭️ Skipping seiffert_bet, already exists: {seiffert_bet}")
    else:
        run(["bet", resampled_seiffert, seiffert_bet, "-f", "0.7"])

#2a) BET sodium TSC datasets which not reference
if radial_tsc_file:
    radial_tsc_bet = f"{strip_ext(radial_tsc_file)}_bet.nii.gz"
    if os.path.exists(radial_tsc_bet):
        print(f"⏭️ Skipping radial_tsc_bet, already exists: {radial_tsc_bet}")
    else:
        run(["bet", radial_tsc_file, radial_tsc_bet, "-f", "0.7"])

if floret_tsc_file:
    floret_tsc_bet = f"{strip_ext(floret_tsc_file)}_bet.nii.gz"
    if os.path.exists(floret_tsc_bet):
        print(f"⏭️ Skipping floret_tsc_bet, already exists: {floret_tsc_bet}")
    else:
        run(["bet", floret_tsc_file, floret_tsc_bet, "-f", "0.7"])

if seiffert_other_tsc:
    seiffert_tsc_bet = f"{strip_ext(resampled_seiffert_tsc)}_bet.nii.gz"
    if os.path.exists(seiffert_tsc_bet):
        print(f"⏭️ Skipping seiffert_tsc_bet, already exists: {seiffert_tsc_bet}")
    else:
        run(["bet", resampled_seiffert_tsc, seiffert_tsc_bet, "-f", "0.7"])

# 3) Flirt BET outputs to sodium ref -
# This is the main key alignment we want. After this, the outputs will be in same space
# as the reference sodium
if radial_file:
    radial_align = f"{strip_ext(radial_bet)}_align12dof.nii.gz"
    radial_mat   = "radial_to_ref.mat"
    if os.path.exists(radial_align):
        print(f"⏭️ Skipping radial_align, already exists: {radial_align}")
    else:
        run([
            "flirt", "-in", radial_bet, "-ref", sodium_ref_file,
            "-out", radial_align, "-omat", radial_mat
        ])
if floret_file:
    floret_align = f"{strip_ext(floret_bet)}_align12dof.nii.gz"
    floret_mat   = "floret_to_ref.mat"
    if os.path.exists(floret_align):
        print(f"⏭️ Skipping floret_align, already exists: {floret_align}")
    else:
        run([
            "flirt", "-in", floret_bet, "-ref", sodium_ref_file,
            "-out", floret_align, "-omat", floret_mat
        ])
if seiffert_other:
    seiffert_align = f"{strip_ext(seiffert_bet)}_align12dof.nii.gz"
    seiffert_mat   = "seiffert_to_ref.mat"
    if os.path.exists(seiffert_align):
        print(f"⏭️ Skipping seiffert_align, already exists: {seiffert_align}")
    else:
        run([
            "flirt", "-in", seiffert_bet, "-ref", sodium_ref_file,
            "-out", seiffert_align, "-omat", seiffert_mat
        ])


#4 ) Flirt non BET outputs to sodium ref
if radial_file:
    radial_align_nobet = f"{strip_ext(radial_file)}_align12dof.nii.gz"
    if os.path.exists(radial_align_nobet):
        print(f"⏭️ Skipping radial_align_nobet, already exists: {radial_align_nobet}")
    else:
        run([
            "flirt", "-in", radial_file, "-ref", sodium_ref_file,
            "-out", radial_align_nobet, "-init", radial_mat, "-applyxfm"
        ])
if floret_file:
    floret_align_nobet = f"{strip_ext(floret_file)}_align12dof.nii.gz"
    if os.path.exists(floret_align_nobet):
        print(f"⏭️ Skipping floret_align_nobet, already exists: {floret_align_nobet}")
    else:
        run([
            "flirt", "-in", floret_file, "-ref", sodium_ref_file,
            "-out", floret_align_nobet, "-init", floret_mat, "-applyxfm"
        ])
if seiffert_other:
    seiffert_align_nobet = f"{strip_ext(resampled_seiffert)}_align12dof.nii.gz"
    if os.path.exists(seiffert_align_nobet):
        print(f"⏭️ Skipping seiffert_align_nobet, already exists: {seiffert_align_nobet}")
    else:
        run([
            "flirt", "-in", resampled_seiffert, "-ref", sodium_ref_file,
            "-out", seiffert_align_nobet, "-init", seiffert_mat, "-applyxfm"
        ])


#5 ) Now we apply transforms to TSC files.

# if radial_tsc_file:
#     print(f"✅ Radial file: {radial_tsc_file}")
# else:
#     print("ℹ️ No radial file (probably reference or missing)")
# if floret_tsc_file:
#     print(f"✅ Floret file: {floret_tsc_file}")
# else:
#     print("ℹ️ No floret file (probably reference or missing)")

if radial_tsc_file:
    radial_tsc_align = f"{strip_ext(radial_tsc_file)}_align12dof.nii.gz"
    if os.path.exists(radial_tsc_align):
        print(f"⏭️ Skipping radial_tsc_align, already exists: {radial_tsc_align}")
    else:
        run([
            "flirt", "-in", radial_tsc_file, "-ref", sodium_ref_file,
            "-out", radial_tsc_align, "-init", radial_mat, "-applyxfm"
        ])
if floret_tsc_file:
    floret_tsc_align = f"{strip_ext(floret_tsc_file)}_align12dof.nii.gz"
    if os.path.exists(floret_tsc_align):
        print(f"⏭️ Skipping floret_tsc_align, already exists: {floret_tsc_align}")
    else:
        run([
            "flirt", "-in", floret_tsc_file, "-ref", sodium_ref_file,
            "-out", floret_tsc_align, "-init", floret_mat, "-applyxfm"
        ])
if seiffert_other_tsc:
    seiffert_tsc_align = f"{strip_ext(resampled_seiffert_tsc)}_align12dof.nii.gz"
    if os.path.exists(seiffert_tsc_align):
        print(f"⏭️ Skipping seiffert_tsc_align, already exists: {seiffert_tsc_align}")
    else:
        run([
            "flirt", "-in", resampled_seiffert_tsc, "-ref", sodium_ref_file,
            "-out", seiffert_tsc_align, "-init", seiffert_mat, "-applyxfm"
        ])

#5a)

if radial_tsc_file:
    radial_tsc_bet_align = f"{strip_ext(radial_tsc_bet)}_align12dof.nii.gz"
    if os.path.exists(radial_tsc_bet_align):
        print(f"⏭️ Skipping radial_tsc_bet_align, already exists: {radial_tsc_bet_align}")
    else:
        run([
            "flirt", "-in", radial_tsc_bet, "-ref", sodium_ref_file,
            "-out", radial_tsc_bet_align, "-init", radial_mat, "-applyxfm"
        ])
if floret_tsc_file:
    floret_tsc_bet_align = f"{strip_ext(floret_tsc_bet)}_align12dof.nii.gz"
    if os.path.exists(floret_tsc_bet_align):
        print(f"⏭️ Skipping floret_tsc_bet_align, already exists: {floret_tsc_bet_align}")
    else:
        run([
            "flirt", "-in", floret_tsc_bet, "-ref", sodium_ref_file,
            "-out", floret_tsc_bet_align, "-init", floret_mat, "-applyxfm"
        ])
if seiffert_other_tsc:
    seiffert_tsc_bet_align = f"{strip_ext(seiffert_tsc_bet)}_align12dof.nii.gz"
    if os.path.exists(seiffert_tsc_bet_align):
        print(f"⏭️ Skipping seiffert_tsc_bet_align, already exists: {seiffert_tsc_bet_align}")
    else:
        run([
            "flirt", "-in", seiffert_tsc_bet, "-ref", sodium_ref_file,
            "-out", seiffert_tsc_bet_align, "-init", seiffert_mat, "-applyxfm"
        ])


### mprage stuff

mprage_optibrain = os.path.join(ARG1, f"{subject}_MPRAGE_optibrain.nii.gz")
mprage_optibrain_mask = os.path.join(ARG1, f"{subject}_MPRAGE_optibrain_mask.nii.gz")
mprage_optibrain_fast = os.path.join(ARG1, f"{subject}_MPRAGE_optibrain_pve_0.nii.gz")

if not os.path.exists(mprage_optibrain):

    run(["sh", OPTIBET_PATH, "-i", mprage_file])

    #base = os.path.splitext(mprage_file)[0]  # remove .nii or .nii.gz
    base = strip_ext(mprage_file) # remove .nii or .nii.gz
    #print(base)

    #sys.exit(0)
    optibet_brain = f"{base}_optiBET_brain.nii.gz"
    optibet_mask = f"{base}_optiBET_brain_mask.nii.gz"

    # Rename/move to desired output names
    shutil.move(optibet_brain, mprage_optibrain)
    shutil.move(optibet_mask, mprage_optibrain_mask)
    print(f"✅ optiBET brain created: {mprage_optibrain}")

else:
    print(f"⏭️ optiBET brain {mprage_optibrain} already exists, skipping.")


if not os.path.exists(mprage_optibrain_fast):
    run([FASTPATH, mprage_optibrain])
else:
    print(f"✅ Already run FSL FAST on optibrain: {mprage_optibrain_fast}")


########### Need to move data now
#####

subject_root = os.path.dirname(os.path.dirname(ARG3))  # go up from .../pipeline/other_sodium → .../site2
output_dir = os.path.join(subject_root, "outputs")
os.makedirs(output_dir, exist_ok=True)
print(f"📦 Collecting outputs in: {output_dir}")

# print(sodium_ref_file)
# print(sodium_ref_tsc_file)

# Gather files if variable exists and is not None
candidates = ["sodium_ref_file", "sodium_ref_tsc_file",
    "radial_file", "radial_tsc_file",
    "floret_file", "floret_tsc_file",
    "resampled_seiffert", "resampled_seiffert_tsc",
    "radial_align", "floret_align", "seiffert_align",
    "radial_align_nobet", "floret_align_nobet", "seiffert_align_nobet",
    "radial_tsc_align", "floret_tsc_align", "seiffert_tsc_align",
    "radial_tsc_bet_align", "floret_tsc_bet_align", "seiffert_tsc_bet_align",
    # MPRAGE outputs
    "mprage_optibrain", "mprage_optibrain_mask"
]

for var in candidates:
    f = locals().get(var)  # check if variable is defined
    if f and os.path.exists(f):
        dest = os.path.join(output_dir, os.path.basename(f))
        if os.path.exists(dest):
            print('skipping')
        else:
            shutil.copy(f, dest)
        print(f"✅ Copied {os.path.basename(f)} → outputs/")
    else:
        print(f"⚠️ Missing or undefined: {var}")

