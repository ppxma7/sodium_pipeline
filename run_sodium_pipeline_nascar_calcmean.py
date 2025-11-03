#!/usr/bin/env python3
import os
import glob
import subprocess
import sys
from collections import defaultdict
import re

# -------- Configuration --------
ROOTDIR = "/Volumes/nemosine/SAN/NASCAR"
subjects = ["Subject1", "Subject2", "Subject3", "Subject4"]
sites = ["site1", "site2", "site3"]

OUTDIR = os.path.join(ROOTDIR, "groupstats")
os.makedirs(OUTDIR, exist_ok=True)
OUTDIR2 = os.path.join(ROOTDIR, "groupstats_acrosssub")
os.makedirs(OUTDIR2, exist_ok=True)

MISSING_LOG = os.path.join(OUTDIR, "missing_files.log")
MISSING_LOG2 = os.path.join(OUTDIR2, "missing_files.log")


# Helper: write missing info to log, choose log file
def log_missing(msg, logfile=MISSING_LOG):
    with open(logfile, "a") as fh:
        fh.write(msg + "\n")



def get_scanbase(filename):
    """
    Take a full filename and return a normalised scanbase for grouping.
    - Ignore BET files
    - Strip align suffixes and final MNI suffix
    """
    base = os.path.basename(filename)

    # Ignore BET files altogether
    if "_bet" in base:
        return None

    # Remove the trailing MNI bit
    #base = re.sub(r'_toMPRAGE_MNI\.nii\.gz$', '', base)
    base = re.sub(r'_alignedtoRef_toMPRAGE_MNI\.nii\.gz$', '', base)


    # Remove align/bet suffixes
    #base = re.sub(r'_bet_align12dof$', '', base)
    #base = re.sub(r'_align12dof$', '', base)
    #base = re.sub(r'_bet$', '', base)

    return base

def fsl_mean_std(in_files, out_prefix):
    """Compute mean and std using FSL (fslmerge + fslmaths), skip if already exists."""
    
    mean_file = f"{out_prefix}_mean.nii.gz"
    std_file = f"{out_prefix}_std.nii.gz"
    
    # Skip if both output files already exist
    if os.path.exists(mean_file) and os.path.exists(std_file):
        print(f"⏭️ Skipping {out_prefix}, mean and std already exist.")
        return

    if len(in_files) == 0:
        print(f"⚠️ No files for {out_prefix}")
        return

    if len(in_files) == 1:
        subprocess.run(["cp", in_files[0], mean_file], check=True)
        subprocess.run(["fslmaths", in_files[0], "-mul", "0", std_file], check=True)
        print(f"✅ {out_prefix}: 1 file (copied as mean, zero std)")
        return

    merged_4d = f"{out_prefix}_merged.nii.gz"
    subprocess.run(["fslmerge", "-t", merged_4d] + in_files, check=True)
    subprocess.run(["fslmaths", merged_4d, "-Tmean", mean_file], check=True)
    subprocess.run(["fslmaths", merged_4d, "-Tstd", std_file], check=True)
    try:
        os.remove(merged_4d)
    except OSError:
        pass

    print(f"✅ {out_prefix}: computed mean/std (n={len(in_files)})")


# -------- Build nested dictionary structure --------
by_subject_site = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

pattern = os.path.join(ROOTDIR, "*", "*", "outputs_mni", "*_MNI.nii.gz")
all_files = glob.glob(pattern)

for f in all_files:
    parts = f.split(os.sep)
    if len(parts) < 5:
        print(f"⚠️ Skipping unexpected path: {f}")
        continue

    subj = parts[-4]
    site = parts[-3]

    scanbase = get_scanbase(f)
    if scanbase is None:
        print(f"🚫 Ignoring BET file: {os.path.basename(f)}")
        continue

    by_subject_site[subj][site][scanbase].append(f)
    print(f"✅ Added: subj={subj}, site={site}, scanbase={scanbase}")

print(f"\n🧾 Total files found: {len(all_files)} (after filtering: {sum(len(files) for subj_dict in by_subject_site.values() for site_dict in subj_dict.values() for files in site_dict.values())})")

if len(all_files) == 0:
    sys.exit("❌ No input files found — check paths!")

# Clear previous missing log
if os.path.exists(MISSING_LOG):
    os.remove(MISSING_LOG)

#sys.exit(0)


# -------- 1. Across SITE mean/std for each subject --------
print("\n📊 Step 1: Across sites per subject")

for subj in subjects:
    print(f"\n👤 Subject: {subj}")

    # collect all scanbases present in any site for this subject
    all_scanbases = set()
    for site in sites:
        site_scanbases = by_subject_site[subj][site].keys()
        all_scanbases.update(site_scanbases)

    if not all_scanbases:
        print(f"  ⚠️ No scanbases found for {subj}, skipping.")
        continue

    for scanbase in sorted(all_scanbases):
        print(f"  ▶ Scanbase: {scanbase}")

        # collect files per site and also track if any site is missing
        in_files = []
        missing_sites = []
        for site in sites:
            site_files = by_subject_site[subj][site].get(scanbase, [])
            if not site_files:
                missing_sites.append(site)
                print(f"    - {site}: MISSING")
            else:
                print(f"    - {site}: {len(site_files)} file(s)")
                for fn in site_files:
                    print(f"       {fn}")
            in_files.extend(site_files)

        # If ANY site missing, skip this scanbase for this subject
        if missing_sites or len(in_files) < len(sites):
            msg = (f"SKIP (across-sites) {subj}/{scanbase}: missing sites "
                   f"{', '.join(missing_sites)} (found {len(in_files)} files)")
            print(f"    ⚠️ {msg}")
            log_missing(msg, MISSING_LOG)
            continue

        # All sites present -> proceed
        print(f"    ✅ All sites present. Total files to merge: {len(in_files)}")
        out_prefix = os.path.join(OUTDIR, f"{subj}_{scanbase}_acrossSites")
        # DEBUG STOP (uncomment to stop here for inspection)
        #sys.exit(0)
        print(f"Subject {subj} scanbase '{scanbase}':")
        for site in sites:
            print(f"  {site}: {by_subject_site[subj][site].get(scanbase, [])}")
        #sys.exit(0)

        fsl_mean_std(in_files, out_prefix)

print(f"Missing/skipped items (if any) written to: {MISSING_LOG}")
#sys.exit(0)


# -------- 2. Across SUBJECT mean/std for each site --------
print("\n📊 Step 2: Across subjects per site")

for site in sites:
    print(f"\n🏠 Site: {site}")

    # collect all scanbases present in any subject for this site
    all_scanbases = set()
    for subj in subjects:
        subj_scanbases = by_subject_site[subj][site].keys()
        all_scanbases.update(subj_scanbases)

    if not all_scanbases:
        print(f"  ⚠️ No scanbases found for {site}, skipping.")
        continue

    for scanbase in sorted(all_scanbases):
        print(f"  ▶ Scanbase: {scanbase}")

        in_files = []
        missing_subjects = []
        for subj in subjects:
            subj_files = by_subject_site[subj][site].get(scanbase, [])
            if not subj_files:
                missing_subjects.append(subj)
                print(f"    - {subj}: MISSING")
            else:
                print(f"    - {subj}: {len(subj_files)} file(s)")
                for fn in subj_files:
                    print(f"       {fn}")
            in_files.extend(subj_files)

        # If ANY subject missing, skip this scanbase for this site
        if missing_subjects or len(in_files) < len(subjects):
            msg = (f"SKIP (across-subjects) {site}/{scanbase}: missing subjects "
                   f"{', '.join(missing_subjects)} (found {len(in_files)} files)")
            print(f"    ⚠️ {msg}")
            log_missing(msg,MISSING_LOG2)
            continue

        # All subjects present -> proceed
        print(f"    ✅ All subjects present. Total files to merge: {len(in_files)}")
        out_prefix = os.path.join(OUTDIR2, f"{site}_{scanbase}_acrossSubjects")
        # DEBUG STOP (uncomment to stop here for inspection)
        #sys.exit(0)
        fsl_mean_std(in_files, out_prefix)

print(f"Missing/skipped items (if any) written to: {MISSING_LOG2}")
print("\n🎉 All done.")

