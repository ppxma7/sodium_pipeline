import os
import glob
import argparse
import subprocess
import nibabel as nib
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import csv

FSLDIR = "/usr/local/fsl"

# ---------- ARGPARSE ----------
parser = argparse.ArgumentParser(description="Run sodium MRI atlas")
parser.add_argument("ARG1", help="Basename of reference sodium (e.g. floret, radial, seiffert)")
parser.add_argument("ARG2", help="Site (1 or 2)")
parser.add_argument("ARG3", help="Subject")

args = parser.parse_args()
ARG1 = args.ARG1
ARG2 = args.ARG2
ARG3 = args.ARG3


def load_atlas_labels(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    labels = {}
    for label in root.findall(".//label"):
        idx = int(label.get("index"))
        name = label.text
        labels[idx] = name
    return labels

def strip_ext(fname):
    if fname.endswith(".nii.gz"):
        return fname[:-7]  # strip ".nii.gz"
    else:
        return os.path.splitext(fname)[0]


def roi_table(sodium_file, atlas_file, out_csv, max_roi=47):
    # Load sodium image and atlas
    sodium_img = nib.load(sodium_file)
    atlas_img = nib.load(atlas_file)

    sodium_data = sodium_img.get_fdata()
    atlas_data = atlas_img.get_fdata().astype(int)

    # Get unique ROI labels (limit to Harvard-Oxford range)
    roi_labels = np.unique(atlas_data)
    roi_labels = roi_labels[(roi_labels >= 0) & (roi_labels <= max_roi)]

    # Load ROI names
    xml_file = os.path.join(FSLDIR, "data/atlases/HarvardOxford-Cortical.xml")
    label_dict = load_atlas_labels(xml_file)

    results = []
    for roi in roi_labels:
        mask = atlas_data == roi
        values = sodium_data[mask]

        if values.size > 0:
            mean_val = np.mean(values)
            std_val = np.std(values)
            median_val = np.median(values)
            roi_name = label_dict.get(roi, f"ROI_{roi}")
            results.append([roi, roi_name, mean_val, std_val, median_val])

    df = pd.DataFrame(results, columns=["ROI", "Name", "Mean", "StdDev", "Median"])
    df.to_csv(out_csv, index=False)
    print(f"✅ ROI stats saved to {out_csv}")

def roi_table_catchexceptions(sodium_file, atlas_file, out_csv, max_roi=47):
    try:
        sodium_img = nib.load(sodium_file)
        atlas_img = nib.load(atlas_file)

        sodium_data = sodium_img.get_fdata()
        atlas_data = atlas_img.get_fdata().astype(int)

        # --- Handle shape mismatch automatically ---
        if sodium_data.shape != atlas_data.shape:
            print(f"⚠️ Shape mismatch: {os.path.basename(sodium_file)} "
                  f"{sodium_data.shape} vs atlas {atlas_data.shape} — resampling sodium to atlas space...")

            resampled_sodium = f"{strip_ext(sodium_file)}_resampled_to_atlas.nii.gz"
            subprocess.run([
                f"{FSLDIR}/bin/flirt",
                "-in", sodium_file,
                "-ref", atlas_file,
                "-out", resampled_sodium,
                "-applyxfm",
                "-usesqform"
            ], check=True)

            sodium_file = resampled_sodium
            sodium_img = nib.load(sodium_file)
            sodium_data = sodium_img.get_fdata()
            print(f"✅ Resampled sodium saved as {os.path.basename(resampled_sodium)}")

        # --- ROI extraction as before ---
        roi_labels = np.unique(atlas_data)
        roi_labels = roi_labels[(roi_labels >= 0) & (roi_labels <= max_roi)]

        xml_file = os.path.join(FSLDIR, "data/atlases/HarvardOxford-Cortical.xml")
        label_dict = load_atlas_labels(xml_file)

        results = []
        for roi in roi_labels:
            mask = atlas_data == roi
            values = sodium_data[mask]
            if values.size > 0:
                mean_val = np.mean(values)
                std_val = np.std(values)
                median_val = np.median(values)
                roi_name = label_dict.get(roi, f"ROI_{roi}")
                results.append([roi, roi_name, mean_val, std_val, median_val])

        if not results:
            print(f"⚠️ No valid ROI data found for {os.path.basename(sodium_file)} — skipping.")
            return

        df = pd.DataFrame(results, columns=["ROI", "Name", "Mean", "StdDev", "Median"])
        df.to_csv(out_csv, index=False)
        print(f"✅ ROI stats saved to {out_csv}")

    except Exception as e:
        print(f"❌ Error processing {os.path.basename(sodium_file)}: {e}")
        return




# --- Paths ---
subject = ARG3   # <-- replace or parse dynamically
base_dir = "/Volumes/nemosine/SAN/NASCAR/"
outputs_mni = os.path.join(base_dir, subject, ARG2, "outputs_mni")
outputs_native = os.path.join(base_dir, subject, ARG2, "outputs")
outputs_pve_native = os.path.join(base_dir, subject, ARG2, "outputs_pve_native")


atlas_mni = os.path.join(FSLDIR, "data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz")
atlas_native = os.path.join(outputs_native, f"{subject}_atlas_in_sodium.nii.gz")

# --- Process MNI-space sodiums ---
#mni_sodiums = glob.glob(os.path.join(outputs_mni, "*MNI.nii.gz"))
# include new naming scheme (alignedtoRef_toMPRAGE_MNI)
mni_sodiums = sorted(glob.glob(os.path.join(outputs_mni, "*_alignedtoRef_toMPRAGE_MNI.nii.gz")))

for f in mni_sodiums:
    out_csv = f"{strip_ext(f)}_ROIstats.csv"
    if os.path.exists(out_csv):
        print('skipping')
    else:
        roi_table(f, atlas_mni, out_csv)

####################################################################################



# #--- Process native-space sodiums ---
# # Determine reference sodium (passed from command line)
# reference = ARG1.lower()

# # Define all possible sodium types
# all_sodiums = ["floret", "radial", "seiffert"]

# # The non-reference ones
# others = [s for s in all_sodiums if s != reference]

# # Start list
# native_sodiums = []

# # --- 1. Add reference sodium files (no align12dof) ---
# for suffix in ["", "_TSC", "_2375", "_TSC_2375"]:
#     matches = glob.glob(os.path.join(outputs_native, f"{reference}{suffix}.nii*"))
#     if matches:
#         native_sodiums.extend(matches)
#     else:
#         print(f"⚠️ Missing reference sodium file: {reference}{suffix}.nii*")

# # --- 2. Add non-reference sodiums (alignedtoRef convention) ---
# for name in others:
#     #for suffix in ["", "_TSC", "_2375", "_TSC_2375", "_TSC_3_bottles", "_TSC_3_bottles_2375"]:
#     #for suffix in ["", "_TSC", "_2375", "_TSC_2375", "_TSC_3_bottles", "_TSC_3_bottles_2375", "_TSC_B1"]:
#     for suffix in [
#         "",
#         "_TSC",
#         "_2375",
#         "_TSC_2375",
#         "_TSC_B1",
#         "_TSC_B1_2375",
#         "_TSC_B1_mode",
#         "_TSC_B1_mode_2375",
#         "_TSC_3_bottles",
#         "_TSC_3_bottles_2375",
#     ]:
#         matches = glob.glob(os.path.join(outputs_native, f"{name}{suffix}_alignedtoRef.nii*"))
#         if matches:
#             native_sodiums.extend(matches)
#         else:
#             print(f"⚠️ Missing aligned sodium file: {name}{suffix}_alignedtoRef.nii*")

# print("🧾 Files to process:")
# for f in native_sodiums:
#     print(f"  - {os.path.basename(f)}")


# for f in native_sodiums:
#     if os.path.exists(f):
#         out_csv = f"{strip_ext(f)}_ROIstats.csv"
#         if os.path.exists(out_csv):
#             print("skipping")
#         else:
#             roi_table_catchexceptions(f, atlas_native, out_csv)
#     else:
#         print(f"⚠️ Missing sodium file: {f}")

# --- Process native-space sodiums ---
reference = ARG1.lower()
all_sodiums = ["floret", "radial", "seiffert"]
native_sodiums = []

# Comprehensive suffix list
suffixes = [
    "",
    "_TSC",
    "_2375",
    "_TSC_2375",
    "_TSC_B1",
    "_TSC_B1_2375",
    "_TSC_B1_mode",
    "_TSC_B1_mode_2375",
    "_TSC_3_bottles",
    "_TSC_3_bottles_2375",
    "_TSC_3_bottles_mode",
    "_TSC_3_bottles_mode_2375"
]

# Unified loop: handle reference + others
for name in all_sodiums:
    for suffix in suffixes:
        # For reference sodium, check both raw and aligned versions
        if name == reference:
            patterns = [
                f"{name}{suffix}.nii*",
                f"{name}{suffix}_alignedtoRef.nii*",
            ]
        else:
            # For non-reference sodiums, only alignedtoRef versions make sense
            patterns = [f"{name}{suffix}_alignedtoRef.nii*"]

        for pattern in patterns:
            matches = glob.glob(os.path.join(outputs_native, pattern))
            if matches:
                native_sodiums.extend(matches)
            else:
                print(f"⚠️ Missing sodium file: {pattern}")

print("🧾 Files to process:")
for f in native_sodiums:
    print(f"  - {os.path.basename(f)}")

# --- Run ROI extraction ---
for f in native_sodiums:
    if os.path.exists(f):
        out_csv = f"{strip_ext(f)}_ROIstats.csv"
        if os.path.exists(out_csv):
            print("skipping")
        else:
            roi_table_catchexceptions(f, atlas_native, out_csv)
    else:
        print(f"⚠️ Missing sodium file: {f}")



# now apply atlas code to pve native space
print("\n--- Processing PVE files with atlas ---")

# 1. Collect all PVE files under outputs_pve_native for this subject/site
pve_files = sorted(
    glob.glob(os.path.join(outputs_pve_native, "*pve*.nii*"))
)

if not pve_files:
    print(f"⚠️ No PVE files found in {outputs_pve_native}")
else:
    print(f"🧾 Found {len(pve_files)} PVE files:")
    for f in pve_files:
        print(f"  - {os.path.basename(f)}")


# 2. Compute global stats
all_results = []
for f in pve_files:
    summary_csv = os.path.join(outputs_pve_native, "PVE_global_summary.csv")
    #out_csv = f"{strip_ext(f)}_ROIstats.csv"  # keep your naming convention

    # if os.path.exists(summary_csv):
    #     print(f"⏭️ Skipping (already processed): {os.path.basename(summary_csv)}")
    #     continue

    if not os.path.exists(f):
        print(f"⚠️ Missing PVE file: {f}")
        continue

    try:
        img = nib.load(f)
        data = img.get_fdata()

        # Clean up: drop NaNs and zeros (optional)
        data = data[np.isfinite(data)]
        data = data[data > 0]

        # Compute stats
        mean_val = np.mean(data)
        std_val = np.std(data)
        median_val = np.median(data)
        q75, q25 = np.percentile(data, [75, 25])
        iqr_val = q75 - q25

        # Create DataFrame (easy to append or merge later)
        df = pd.DataFrame({
            "Filename": [os.path.basename(f)],
            "Mean": [mean_val],
            "Std": [std_val],
            "Median": [median_val],
            "IQR": [iqr_val]
        })

        # Save single CSV per file
        all_results.append(df)
        #df.to_csv(out_csv, index=False)
        print(f"✅ Saved global stats → {os.path.basename(out_csv)}")

    except Exception as e:
        print(f"❌ Failed on {f}: {e}")

# Save combined summary at the end
if all_results:
    pd.concat(all_results, ignore_index=True).to_csv(summary_csv, index=False)
    print(f"📊 Saved combined summary → {summary_csv}")
else:
    print("⚠️ No results to save.")
