#!/usr/bin/env python3
"""
BIDS App: Precision Functional Mapping via Template Matching (ReproTM)

Pipeline steps:
  1. Z-score dconn files
  2. Generate precision functional maps via template matching
  3. Clean up networks below minimum size
  4. Convert dscalar to dlabel

Usage:
  python run.py <bids_dir> <output_dir> participant [options]
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

global zscore_script, repro_tm_script, minsize_script
zscore_script=Path('./zscore_dconn/zscore_dconn_v1.0.0.py')
repro_tm_script=Path('./ReproTM/ReproTM_v1.0.0.py')
minsize_script=Path('./minsize/minsize_v1.0.0.py')

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="BIDS App – Precision Functional Mapping via Template Matching (ReproTM)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- BIDS required positional args ---
    p.add_argument("bids_dir", type=Path,
                   help="Root BIDS dataset directory (input dconns live here).")
    p.add_argument("output_dir", type=Path,
                   help="Output directory (derivatives).")
    p.add_argument("analysis_level", choices=["participant"],
                   help="Processing level. Only 'participant' is supported.")

    # --- BIDS participant filtering ---
    p.add_argument("--participant_label", nargs="+", metavar="SUB",
                   help="One or more participant labels to process (without 'sub-' prefix).")
    p.add_argument("--session_label", nargs="+", metavar="SES",
                   help="One or more session labels to process (without 'ses-' prefix).")
    p.add_argument("--task_label", nargs="+", metavar="TASK",
                   help="One or more task labels to process.")

    # --- Pipeline step toggles ---
    steps = p.add_argument_group("Pipeline steps (all enabled by default)")
    steps.add_argument("--skip_zscore", action="store_true",
                       help="Skip step 1: z-scoring dconns.")
    steps.add_argument("--skip_repro_tm", action="store_true",
                       help="Skip step 2: template matching.")
    steps.add_argument("--skip_minsize", action="store_true",
                       help="Skip step 3: minimum-size cleanup.")
    steps.add_argument("--skip_dlabel", action="store_true",
                       help="Skip step 4: dscalar -> dlabel conversion.")

    # --- Code paths ---
    code = p.add_argument_group("Code / tool paths")
    # code.add_argument("--zscore_script", type=Path,
    #                   default=Path("/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/code/zscore_dconn/zscore_dconn_v1.0.0.py"),
    #                   help="Path to zscore_dconn_v1.0.0.py.")
    # code.add_argument("--repro_tm_script", type=Path,
    #                   default=Path("/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/code/ReproTM/ReproTM_v1.0.0.py"),
    #                   help="Path to ReproTM_v1.0.0.py.")
    # code.add_argument("--minsize_script", type=Path,
    #                   default=Path("/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/code/minsize/minsize_v1.0.0.py"),
    #                   help="Path to minsize_v1.0.0.py.")
    code.add_argument("--wb_command", type=Path, required=True,
                      default=Path("wb_command"),
                      help="Path to Connectome Workbench wb_command executable.")

    # --- Template matching options ---
    tm = p.add_argument_group("Template matching (step 2)")
    tm.add_argument("--template_infile", type=Path, required=True,
                    help="Path to the template .mat file.")
    tm.add_argument("--template_networks", type=str, required=True,
                    default="DMN Vis FP NaN DAN NaN VAN Sal AMN SMd SMl Aud Tpole MTL PMN PON NaN SCAN",
                    help="Space-separated list of network names matching the template.")
    tm.add_argument("--surface_only", action="store_true",
                    help="Run template matching for cortex only.")
    tm.add_argument("--template_thresholding", action="store_true",
                    help="Threshold template seedmaps to top connections.")
    tm.add_argument("--template_minthreshold", type=int, default=1,
                    help="Minimum threshold for template thresholding.")
    tm.add_argument("--refineSCAN", action="store_true",
                    help="Re-run template matching with higher threshold for SCAN/SMl/SMd.")
    tm.add_argument("--refineSCAN_minthreshold", type=int, default=3,
                    help="Minimum threshold applied to SCAN, SMl, and SMd networks.")

    # --- Min-size cleanup options ---
    ms = p.add_argument_group("Min-size cleanup (step 3)")
    ms.add_argument("--minsize", type=int, default=30,
                    help="Minimum network size (vertices). Networks smaller than this are removed.")

    # --- Label file for step 4 ---
    lbl = p.add_argument_group("dscalar -> dlabel conversion (step 4)")
    lbl.add_argument("--label_file", type=Path, required=False, default=Path("./dscalar2dlabel/support_files/tpl-ABCC2026-a3-9to16_15networks_label_list.txt"),
                     help="Path to label list .txt file for wb_command -cifti-label-import.")

    # --- Misc ---
    p.add_argument("--dry_run", action="store_true",
                   help="Print commands without executing them.")
    p.add_argument("--verbose", action="store_true",
                   help="Enable verbose/debug logging.")

    return p


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd: list, dry_run: bool = False):
    """Run a shell command, raising on failure."""
    cmd_str = " ".join(str(c) for c in cmd)
    log.info("CMD: %s", cmd_str)
    if dry_run:
        log.info("[dry-run] skipping execution")
        return
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        log.error("Command failed (exit %d): %s", result.returncode, cmd_str)
        sys.exit(result.returncode)


def ensure_dir(path: Path, dry_run: bool = False):
    if not dry_run:
        path.mkdir(parents=True, exist_ok=True)


def discover_subjects(bids_dir: Path, participant_labels) -> list:
    """Return sorted list of subject IDs present in bids_dir."""
    all_subs = sorted(
        d.name.replace("sub-", "") for d in bids_dir.iterdir()
        if d.is_dir() and d.name.startswith("sub-")
    )
    if participant_labels:
        filtered = [s for s in all_subs if s in participant_labels]
        missing = set(participant_labels) - set(all_subs)
        if missing:
            log.warning("Requested participants not found in BIDS dir: %s", missing)
        return filtered
    log.info("No participant labels specified; processing all subjects found in BIDS dir.")
    return all_subs


def discover_sessions(sub_dir: Path, session_labels) -> list:
    """Return sorted list of session IDs for a subject directory."""
    all_ses = sorted(
        d.name.replace("ses-", "") for d in sub_dir.iterdir()
        if d.is_dir() and d.name.startswith("ses-")
    )
    if session_labels:
        filtered = [s for s in all_ses if s in session_labels]
        missing = set(session_labels) - set(all_ses)
        if missing:
            log.warning("Requested sessions not found in subject directory: %s", missing)
        return filtered
    log.info("No session labels specified; processing all sessions found for subject.")
    return all_ses


def discover_tasks(func_dir: Path, task_labels) -> list:
    """Return all dconn.nii files, optionally filtered by task labels."""
    dconn_files = sorted(func_dir.glob("*.dconn.nii"))
    
    if not task_labels:
        log.info("No task labels specified; processing all dconn files found in func directory.")
        # Return all dconn files if no task filter provided
        return dconn_files
    
    # Filter dconn files by task labels
    filtered_files = []
    for f in dconn_files:
        parts = {p.split("-")[0]: p.split("-")[1] for p in f.stem.split("_") if "-" in p}
        if "task" in parts and parts["task"] in task_labels:
            filtered_files.append(f)
    return filtered_files


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def step1_zscore(args, dconn_in, sub, ses, task, bids_dir, output_dir):
    """Z-score a dconn file."""
    log.info("######## Running z-score ########")
    dconn_out_dir = output_dir / f"sub-{sub}" / f"ses-{ses}" / "func"
    dconn_out = dconn_out_dir / f"sub-{sub}_ses-{ses}_task-{task}_stat-zscored.dconn.nii"

    ensure_dir(dconn_out_dir, args.dry_run)

    if not args.dry_run and not dconn_in.exists():
        log.warning("Input dconn not found, skipping: %s", dconn_in)
        return None

    cmd = [
        sys.executable, "-u", str(zscore_script),
        "--dconn_infile", str(dconn_in),
        "--dconn_outfile", str(dconn_out),
    ]
    run_cmd(cmd, args.dry_run)
    return dconn_out


def step2_repro_tm(args, sub, ses, task, template, zscored_dconn, output_dir):
    """Run template matching to generate precision functional maps."""
    out_dir = output_dir / f"sub-{sub}" / f"ses-{ses}" / "func"
    ensure_dir(out_dir, args.dry_run)
    mat_out       = out_dir / f"sub-{sub}_ses-{ses}_task-{task}_ReproTM_template-{template}.mat"
    dscalar_out   = out_dir / f"sub-{sub}_ses-{ses}_task-{task}_ReproTM_template-{template}.dscalar.nii"
    dscalar_scan  = out_dir / f"sub-{sub}_ses-{ses}_task-{task}_ReproTM_template-{template}_refine-SCAN.dscalar.nii"
    log.info("######## Running template matching ########")

    cmd = [
        sys.executable, "-u", str(repro_tm_script),
        "--dconn_infile",    str(zscored_dconn),
        "--template_infile", str(args.template_infile),
        "--template_networks", *args.template_networks.split(),
        "--mat_outfile",          str(mat_out),
        "--dscalar_outfile",      str(dscalar_out),
        "--dscalarSCANrefined_outfile", str(dscalar_scan),
    ]

    if args.surface_only:
        cmd.append("--surface_only")
    if args.template_thresholding:
        cmd += ["--template_thresholding",
                "--template_minthreshold", str(args.template_minthreshold)]
    if args.refineSCAN:
        cmd += ["--refineSCAN",
                "--refineSCAN_minthreshold", str(args.refineSCAN_minthreshold)]

    run_cmd(cmd, args.dry_run)
    return dscalar_out, dscalar_scan


def step3_minsize(args, sub, ses, dscalar_in, output_dir):
    """Remove networks below minimum vertex size."""
    log.info("######## Running Python min-size cleanup ########")
    out_dir = output_dir / f"sub-{sub}" / f"ses-{ses}" / "func"
    ensure_dir(out_dir, args.dry_run)

    filename_stem = dscalar_in.stem.replace(".dscalar", "")

    dscalar_out = out_dir / \
        f"{filename_stem}_minsize-{args.minsize}.dscalar.nii"

    cmd = [
        sys.executable, str(minsize_script),
        "--dscalar_infile",  str(dscalar_in),
        "--dscalar_outfile", str(dscalar_out),
        "--minsize",         str(args.minsize),
    ]
    log.info("############ Cleaning networks below %d threshold ############", args.minsize)
    run_cmd(cmd, args.dry_run)
    return dscalar_out

def step3_minsize_matlab(args, sub, ses, dscalar_in, output_dir):
    """Remove networks below minimum vertex size using MATLAB."""
    log.info("######## Running MATLAB minsize cleanup ########")
    out_dir = output_dir / f"sub-{sub}" / f"ses-{ses}" / "func"
    ensure_dir(out_dir, args.dry_run)

    filename_stem = dscalar_in.stem.replace(".dscalar", "")
    dscalar_out = out_dir / f"{filename_stem}_recolored_minsize{args.minsize}.dscalar.nii"

    log.info("############ Cleaning networks below %d threshold ############", args.minsize)

    # Build MATLAB command
    matlab_addpath = "addpath(genpath('/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/code/minsize'))"
    matlab_func_call = f"clean_dscalars_by_size_simple_v2('{str(dscalar_in)}',[],[],[],[],{args.minsize},[],0,0,1);"
    matlab_command = f"{matlab_addpath};{matlab_func_call}"
    
    matlab_exec = Path("/common/software/install/migrated/matlab/R2019a/bin/matlab")

    cmd = [
        str(matlab_exec),
        "-nodisplay",
        "-nosplash",
        "-batch",
        matlab_command,
    ]
    # log.info("MATLAB minsize command: %s", matlab_command)
    run_cmd(cmd, args.dry_run)
    
    return dscalar_out

def step4_dlabel(args, sub, ses, dscalar_in, output_dir):
    """Convert dscalar to dlabel using wb_command."""
    out_dir = output_dir / f"sub-{sub}" / f"ses-{ses}" / "func"
    ensure_dir(out_dir, args.dry_run)

    stem = dscalar_in.name.replace(".dscalar.nii", "")
    dlabel_out = out_dir / f"{stem}.dlabel.nii"

    cmd = [
        str(args.wb_command),
        "-cifti-label-import",
        str(dscalar_in),
        str(args.label_file),
        str(dlabel_out),
    ]
    run_cmd(cmd, args.dry_run)
    return dlabel_out


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_args(args):
    errors = []
    
    if not args.skip_zscore and not zscore_script.exists():
        errors.append(f"zscore_script not found: {zscore_script}")

    if not args.skip_repro_tm:
        if not repro_tm_script.exists():
            errors.append(f"repro_tm_script not found: {repro_tm_script}")
        if args.template_infile is None:
            errors.append("--template_infile is required for step 2 (template matching).")
        elif not args.template_infile.exists():
            errors.append(f"template_infile not found: {args.template_infile}")

    if not args.skip_minsize and not minsize_script.exists():
        errors.append(f"minsize_script not found: {minsize_script}")

    if not args.skip_dlabel:
        if args.label_file is None:
            errors.append("--label_file is required for step 4 (dscalar -> dlabel).")
        elif not args.label_file.exists():
            errors.append(f"label_file not found: {args.label_file}")

    if errors:
        for e in errors:
            log.error(e)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.dry_run:
        log.info("=== DRY RUN MODE – no commands will be executed ===")

    validate_args(args)

    bids_dir: Path = args.bids_dir.resolve()
    output_dir: Path = args.output_dir.resolve()

    if not bids_dir.exists():
        log.error("BIDS directory not found: %s", bids_dir)
        sys.exit(1)

    ensure_dir(output_dir, args.dry_run)

    subjects = discover_subjects(bids_dir, args.participant_label)
    if not subjects:
        log.error("No matching subjects found in %s", bids_dir)
        sys.exit(1)

    log.info("Subjects to process: %s", subjects)

    for sub in subjects:
        sub_dir = bids_dir / f"sub-{sub}"
        sessions = discover_sessions(sub_dir, args.session_label)

        if not sessions:
            log.warning("No sessions found for sub-%s – skipping.", sub)
            continue

        for ses in sessions:
            func_dir = sub_dir / f"ses-{ses}" / "func"
            if not func_dir.exists():
                log.warning("No func dir for sub-%s ses-%s – skipping.", sub, ses)
                continue

            dconn_files = discover_tasks(func_dir, args.task_label)
            if not dconn_files:
                log.warning("No dconn files for sub-%s ses-%s – skipping.", sub, ses)
                continue

            for dconn_file in dconn_files:
                # Extract task name from filename
                parts = {p.split("-")[0]: p.split("-")[1] for p in dconn_file.stem.split("_") if "-" in p}
                task = parts.get("task", "unknown")
                log.info("--- Processing sub-%s ses-%s task-%s ---", sub, ses, task)

                # ── Step 1: z-score ──────────────────────────────────────
                if not args.skip_zscore:
                    zscored = step1_zscore(args, dconn_file, sub, ses, task, bids_dir, output_dir)
                    if zscored is None:
                        continue
                else:
                    zscored = (output_dir / f"sub-{sub}" / f"ses-{ses}" / "func" /
                               f"sub-{sub}_ses-{ses}_task-{task}_stat-zscored.dconn.nii")
                    log.info("Skipping step 1; expecting z-scored dconn at: %s", zscored)
                
                if "tpl-" in args.template_infile.name:
                    template = args.template_infile.name.split('.')[0].split('tpl-')[-1].split('_')[0]  # get template name without extension for output filenames
                elif "seedmaps-" in args.template_infile.name:
                    template = args.template_infile.name.split('.')[0].split('seedmaps-')[-1].split('_')[0]
                else:
                    template = args.template_infile.stem  # fallback to full stem if no recognizable pattern


                # ── Step 2: template matching ─────────────────────────────
                if not args.skip_repro_tm:
                    dscalar, dscalar_scan = step2_repro_tm(
                        args, sub, ses, task, template, zscored, output_dir)
                else:
                    dscalar = (output_dir / f"sub-{sub}" / f"ses-{ses}" / "func" /
                               f"sub-{sub}_ses-{ses}_task-{task}_ReproTM_template-{template}.dscalar.nii")
                    if args.refineSCAN:
                        dscalar_scan = (output_dir / f"sub-{sub}" / f"ses-{ses}" / "func" /
                                        f"sub-{sub}_ses-{ses}_task-{task}_ReproTM_template-{template}_refine-SCAN.dscalar.nii")
                    log.info("Skipping step 2; expecting dscalar at: %s", dscalar)

                # ── Step 3: min-size cleanup ──────────────────────────────
                if not args.skip_minsize:
                    dscalar_clean = step3_minsize(
                        args, sub, ses, dscalar, output_dir)
                    if args.refineSCAN:
                        log.info("RefineSCAN enabled; using SCAN-refined dscalar for min-size cleanup.")
                        dscalar_clean_scan = step3_minsize(
                            args, sub, ses, dscalar_scan, output_dir)
                else:
                    dscalar_clean = dscalar
                    log.info("Skipping step 3 (min-size cleanup).")

                # ── Step 4: dscalar -> dlabel ─────────────────────────────
                if not args.skip_dlabel:
                    step4_dlabel(args, sub, ses, dscalar_clean, output_dir)
                    if args.refineSCAN:
                        step4_dlabel(args, sub, ses, dscalar_clean_scan, output_dir)
                else:
                    log.info("Skipping step 4 (dlabel conversion).")

    log.info("Pipeline complete.")


if __name__ == "__main__":
    main()