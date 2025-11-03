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
import tempfile


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

FORCE_DIRECT_FLIRT = False  # set True to run direct FLIRT on non-BET images


def get_extension(fname):
    """Return file extension (.nii.gz or .nii) even if input is a list."""
    if isinstance(fname, list):
        if len(fname) == 0:
            return None
        fname = fname[0]  # take the first element

    if fname.endswith(".nii.gz"):
        return ".nii.gz"
    return os.path.splitext(fname)[1]

def run(cmd, check=True):
    print("🔧 Running:", " ".join(cmd))
    subprocess.run(cmd, check=check)

def strip_ext(fname):
    if fname.endswith(".nii.gz"):
        return fname[:-7]  # strip 7 chars for '.nii.gz'
    else:
        return os.path.splitext(fname)[0]
        

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
print(f"This is the sodium reference: {sodium_ref_file}")

# 3. Reference sodium TSC
ref_base = os.path.splitext(os.path.basename(sodium_ref_file))[0]
sodium_ref_tsc_matches = glob.glob(os.path.join(ARG2, f"{ref_base}_TSC*.nii*"))

# sodium_ref_tsc_matches = (
#     glob.glob(os.path.join(ARG2, f"{ref_base}_TSC.nii*")) +
#     glob.glob(os.path.join(ARG2, f"{ref_base}_TSC_3_bottles.nii*"))
# )
if not sodium_ref_tsc_matches:
    print(f"No sodium ref TSC files found in {ARG2} for {ref_base}, continuing without TSC")
    sodium_ref_tsc_files = []
else:
    print(f"Found sodium ref TSC files: {sodium_ref_tsc_matches}")
    sodium_ref_tsc_files = sodium_ref_tsc_matches

# --- Self-register reference sodium (for consistent headers) ---
print(f"🔁 Self-registering reference sodium to itself for header consistency")

def self_register_to_ref(ref_file, tsc_files):
    all_files = [ref_file] + tsc_files
    out_files = []
    for src in all_files:
        fname = os.path.basename(src)
        # ✅ Skip already aligned files
        if "_alignedtoRef" in fname:
            print(f"⏭️ Skipping (already aligned): {fname}")
            continue

        base = strip_ext(src)
        out_file = f"{base}_alignedtoRef.nii.gz"
        if os.path.exists(out_file):
            print(f"⏭️ Skipping (already exists): {fname}")
            continue

        print(f"🟢 Self-registering {fname} to itself")
        run([
            "flirt",
            "-in", src,
            "-ref", ref_file,
            "-out", out_file,
            "-applyxfm",
            "-usesqform"
        ])
        out_files.append(out_file)
    return out_files


# Apply to reference sodium and its TSCs
ref_aligned = self_register_to_ref(sodium_ref_file, sodium_ref_tsc_files)

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


# --- Find TSC files for the two "other sodium" images ---
# other_tsc_files = {}
# for name, f in valid_other_sodiums.items():
#     base = f[:-7] if f.endswith(".nii.gz") else os.path.splitext(f)[0]
#     tsc_matches = glob.glob(f"{base}_TSC.nii*")
#     if len(tsc_matches) == 0:
#         print(f"No TSC file found for {name} ({f}), continuing without TSC")
#         other_tsc_files[name] = None
#     elif len(tsc_matches) > 1:
#         print(f"Multiple TSC files for {name}, using first: {tsc_matches[0]}")
#         other_tsc_files[name] = tsc_matches[0]
#     else:
#         other_tsc_files[name] = tsc_matches[0]

# --- Find TSC files for the "other sodium" images ---
other_tsc_files = {}
for name, f in valid_other_sodiums.items():
    base = f[:-7] if f.endswith(".nii.gz") else os.path.splitext(f)[0]
    tsc_matches = sorted(glob.glob(f"{base}_TSC*.nii*"))

    if not tsc_matches:
        print(f"No TSC files found for {name} ({f}), continuing without TSC")
        other_tsc_files[name] = []
    else:
        print(f"Found TSC files for {name}: {tsc_matches}")
        other_tsc_files[name] = tsc_matches

print("\nReference TSC files:", sodium_ref_tsc_files)
print("Other sodium TSC files:", other_tsc_files)




#sys.exit(0)

######################
# Okay here we have all our files. the sodium images (file1 and file2)
# The reference sodium image
# the 3dt1

######################


# Step 1) Resample Seiffert if present (reference or other)
resampled_seiffert = None
resampled_seiffert_tsc = []


def resample_all_seiffert_tscs(folder):
    """Find and resample all Seiffert TSC files (_TSC*, including 3_bottles)."""
    tscs = sorted(glob.glob(os.path.join(folder, "seiffert_TSC*.nii*")))
    tscs = [t for t in tscs if "_2375" not in os.path.basename(t)]
    resampled = []
    for tsc in tscs:
        resampled_tsc = f"{strip_ext(tsc)}_2375.nii.gz"
        if os.path.exists(resampled_tsc):
            print(f"⏭️ Skipping Seiffert TSC resample, already exists: {resampled_tsc}")
        else:
            run(["flirt", "-in", tsc, "-ref", tsc,
                 "-out", resampled_tsc, "-applyisoxfm", "2.375"])
        resampled.append(resampled_tsc)
    return resampled


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

    # Handle multiple TSCs
    # Resample ALL Seiffert TSCs in the same folder
    ref_dir = os.path.dirname(sodium_ref_file)

    #print(ref_dir)
    #sys.exit(0)
    resampled_seiffert_tscs = resample_all_seiffert_tscs(ref_dir)

    # repoint the reference
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

    # Resample ALL Seiffert TSCs in the same folder
    other_dir = os.path.dirname(seiffert_file)
    resampled_seiffert_tscs = resample_all_seiffert_tscs(other_dir)

else:
    print("⏭️ No Seiffert sodium file found (neither reference nor other), skipping resampling.")




# BET input files
radial_file = other_sodium_files.get("radial")  # single file, OK
floret_file = other_sodium_files.get("floret")  # single file, OK
seiffert_other = other_sodium_files.get("seiffert")  # list of Seiffert files?

# TSC files (use resampled lists for Seiffert)
radial_tsc_file = other_tsc_files.get("radial")  # single
floret_tsc_file = other_tsc_files.get("floret")  # single
# use resampled_seiffert_tscs list from previous step
seiffert_other_tscs = resampled_seiffert_tscs  # list of all resampled Seiffert TSCs

# print(radial_file)
# print(floret_file)
# print(seiffert_other)

# print(radial_tsc_file)
# print(floret_tsc_file)
# print(seiffert_other_tscs)

# sys.exit(0)

def run_bet(files):
    """Run BET on a list of files, skipping any that are already BET or aligned."""
    if isinstance(files, str):
        files = [files]

    bet_files = []
    for f in files:
        base = strip_ext(f)
        fname = os.path.basename(f)

        # Skip already processed files
        if "_bet" in fname or "_align12dof" in fname or "_alignedtoRef" in fname:
            print(f"ℹ️ Skipping BET (already BETed or aligned): {f}")
            bet_files.append(f)
            continue

        out_file = f"{base}_bet.nii.gz"
        if os.path.exists(out_file):
            print(f"⏭️ Skipping BET, already exists: {out_file}")
        else:
            print(f"🔧 Running BET: {f} -> {out_file}")
            run(["bet", f, out_file, "-f", "0.7"])

        bet_files.append(out_file)

    return bet_files



# --- FLIRT section ---
def run_flirt(files, reference_file, apply_xfm=False):
    """
    Run FSL FLIRT alignment to reference.
    - If apply_xfm=True, expects a corresponding *_to_ref.mat file and applies transform.
    - Otherwise, runs full 12-dof alignment.
    Returns list of aligned filenames.
    """
    if isinstance(files, str):
        files = [files]

    aligned_files = []
    for f in files:
        fname = os.path.basename(f)

        # ✅ Skip if already aligned (prevents repeated re-alignment)
        if "_alignedtoRef" in fname:
            print(f"ℹ️ Skipping {f}, looks like a FLIRT output already")
            continue

        base = strip_ext(f)
        align_file = f"{base}_alignedtoRef.nii.gz"      # ✅ new consistent naming
        mat_file   = f"{base}_to_ref.mat"

        if os.path.exists(align_file):
            print(f"⏭️ Skipping FLIRT, already exists: {align_file}")
        else:
            if apply_xfm and os.path.exists(mat_file):
                print(f"🟢 Applying existing transform: {f} -> {align_file}")
                run([
                    "flirt",
                    "-in", f,
                    "-ref", reference_file,
                    "-out", align_file,
                    "-init", mat_file,
                    "-applyxfm"
                ])
            else:
                print(f"🟢 Running full FLIRT: {f} -> {align_file}")
                run([
                    "flirt",
                    "-in", f,
                    "-ref", reference_file,
                    "-out", align_file,
                    "-omat", mat_file
                ])

        aligned_files.append(align_file)

    return aligned_files


# Ensure TSC variables are lists of strings
radial_tsc_list   = radial_tsc_file if isinstance(radial_tsc_file, list) else [radial_tsc_file] if radial_tsc_file else []
floret_tsc_list   = floret_tsc_file if isinstance(floret_tsc_file, list) else [floret_tsc_file] if floret_tsc_file else []
seiffert_tsc_list = seiffert_other_tscs if seiffert_other_tscs else []

# --- BET section ---
radial_bet_files   = run_bet([radial_file]) if radial_file else []
floret_bet_files   = run_bet([floret_file]) if floret_file else []
seiffert_bet_files = run_bet(seiffert_tsc_list if seiffert_tsc_list else [resampled_seiffert]) if resampled_seiffert or seiffert_tsc_list else []

# BET TSCs (now safely flattened)
radial_tsc_bet_files   = run_bet(radial_tsc_list) if radial_tsc_list else []
floret_tsc_bet_files   = run_bet(floret_tsc_list) if floret_tsc_list else []
seiffert_tsc_bet_files = run_bet(seiffert_tsc_list) if seiffert_tsc_list else []

# Step 2a: FLIRT 12-dof on BET images (estimate transforms)
# radial_bet_align   = run_flirt(radial_bet_files, sodium_ref_file)
# floret_bet_align   = run_flirt(floret_bet_files, sodium_ref_file)
# seiffert_bet_align = run_flirt(seiffert_bet_files, sodium_ref_file)

# Step 2a: Estimate transform using BET images (save .mat only)
for f in radial_bet_files + floret_bet_files + seiffert_bet_files:
    mat_file = f"{strip_ext(f)}_to_ref.mat"
    if not os.path.exists(mat_file):
        tmp_out = f"{strip_ext(f)}_tmpflirt.nii.gz"
        print(f"🧠 Estimating transform only (BET): {f}")
        run([
            "flirt",
            "-in", f,
            "-ref", sodium_ref_file,
            "-omat", mat_file,
            "-out", tmp_out
        ])
        os.remove(tmp_out)
    else:
        print(f"⏭️ Transform already exists for: {f}")


# Step 2b: Apply transforms to original + TSC files (using *_bet_to_ref.mat)
def apply_transform(src_files, ref_file, bet_file_list):
    """Apply transform estimated from BET file(s) to original/TSC files."""
    if not bet_file_list:
        print("⚠️ No BET file provided, skipping transform application.")
        return []
    mat_file = strip_ext(bet_file_list[0]) + "_to_ref.mat"
    if not os.path.exists(mat_file):
        print(f"⚠️ Missing transform: {mat_file}")
        return []

    if isinstance(src_files, str):
        src_files = [src_files]

    out_files = []
    for src in src_files:
        fname = os.path.basename(src)
        # Skip BETed or already aligned files
        if any(tag in fname for tag in ["_bet", "_alignedtoRef", "_align12dof"]):
            print(f"ℹ️ Skipping BETed/aligned file: {src}")
            continue

        base = strip_ext(src)
        out_file = f"{base}_alignedtoRef.nii.gz"
        if os.path.exists(out_file):
            print(f"⏭️ Skipping (already aligned): {out_file}")
        else:
            print(f"🟢 Applying transform: {src} using {mat_file}")
            run([
                "flirt",
                "-in", src,
                "-ref", ref_file,
                "-out", out_file,
                "-init", mat_file,
                "-applyxfm"
            ])
        out_files.append(out_file)
    return out_files




# --- Optional: Direct robust FLIRT (skip BET) ---
if FORCE_DIRECT_FLIRT:
    print("⚡ FORCE_DIRECT_FLIRT is ON — running direct robust FLIRT on non-BET images")

    def run_direct_flirt(files, ref_file):
        """Run robust FLIRT directly on full-head inputs (no BET dependency)."""
        if isinstance(files, str):
            files = [files]

        for f in files:
            fname = os.path.basename(f)
            if any(tag in fname for tag in ["_alignedtoRef", "_bet", "_align12dof"]):
                print(f"⏭️ Skipping already aligned/BET file: {fname}")
                continue

            out_file = f"{strip_ext(f)}_alignedtoRef.nii.gz"
            if os.path.exists(out_file):
                print(f"⏭️ Skipping existing: {out_file}")
                continue

            mat_file = f"{strip_ext(f)}_to_ref.mat"
            print(f"🧭 Direct FLIRT: {f} → {out_file}")
            run([
                "flirt",
                "-in", f,
                "-ref", ref_file,
                "-out", out_file,
                "-omat", mat_file,
                "-bins", "256",
                "-cost", "corratio",
                "-searchrx", "-90", "90",
                "-searchry", "-90", "90",
                "-searchrz", "-90", "90",
                "-dof", "12",
                "-interp", "trilinear"
            ])

    # Example: only rerun Seiffert
    run_direct_flirt([resampled_seiffert] + seiffert_tsc_list, sodium_ref_file)
else:
    # Apply to originals and TSCs
    radial_apply   = apply_transform([radial_file] + radial_tsc_list, sodium_ref_file, radial_bet_files)
    floret_apply   = apply_transform([floret_file] + floret_tsc_list, sodium_ref_file, floret_bet_files)
    seiffert_apply = apply_transform([resampled_seiffert] + seiffert_tsc_list, sodium_ref_file, seiffert_bet_files)




#sys.exit(0)
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


########### Collect aligned + MPRAGE outputs ###########

subject_root = os.path.dirname(os.path.dirname(ARG3))  # e.g. .../Subject1/site2
output_dir = os.path.join(subject_root, "outputs")
os.makedirs(output_dir, exist_ok=True)
print(f"📦 Collecting outputs in: {output_dir}")

# 1️⃣ Find all aligned files in the current pipeline tree
aligned_files = []
for root, _, files in os.walk(subject_root):
    for fname in files:
        if fname.endswith("_alignedtoRef.nii.gz"):
            aligned_files.append(os.path.join(root, fname))

# 2️⃣ Explicitly add the MPRAGE outputs
mprage_files = []
for var in ["mprage_optibrain", "mprage_optibrain_mask"]:
    f = locals().get(var)
    if f and os.path.exists(f):
        mprage_files.append(f)

# 3️⃣ Combine and copy to outputs folder
all_to_copy = aligned_files + mprage_files

if not all_to_copy:
    print("⚠️ No aligned or MPRAGE files found to copy.")
else:
    for f in all_to_copy:
        dest = os.path.join(output_dir, os.path.basename(f))
        if os.path.exists(dest):
            print(f"⏭️ Skipping (already exists): {os.path.basename(f)}")
        else:
            shutil.copy(f, dest)
            print(f"✅ Copied {os.path.basename(f)} → outputs/")


