# ReproTM: Reproducible Template Matching

ReproTM (Godfrey et al., in prep) is a set of python tools to facilitate reproducible generation of individualized, person-specific, functional network maps using the template matching network detection algorithm ([Gordon et al., 2017](https://www.sciencedirect.com/science/article/pii/S089662731730613X))

## Overview

This pipeline takes CIFTI dense connectome files (`.dconn.nii`) organized in BIDS format and produces functional network parcellation maps for each participant. It runs four sequential steps:

1. **Z-score** dconn files
2. **Template matching** (ReproTM) to generate individual functional network maps
3. **Minimum size cleanup** clean up network clusters below a greyordinate threshold
4. **Convert** dscalar to dlabel via `wb_command`

Each step can be skipped independently.

>[!CAUTION]
>Running ReproTM as a BIDSapp with run.py is still under development, current recommendation is to manually run each ReproTM step:
>1. */code/zscore_dconn/zscore_dconn.sh
>2. */code/ReproTM/ReproTM.sh
>3. */code/minsize/minsize.sh

## Requirements

- Python 3
- [Nibabel](https://nipy.org/nibabel/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench) (`wb_command`)
- The following scripts present relative to `run.py`:
  - `zscore_dconn/zscore_dconn_v1.0.0.py`
  - `ReproTM/ReproTM_v1.0.0.py`
  - `minsize/minsize_v1.0.0.py`


## Input Data

The pipeline expects a BIDS-formatted input directory with the following structure:

```
bids_dir/
└── sub-<label>/
    └── ses-<label>/
        └── func/
            └── sub-<label>_ses-<label>_task-<label>_*.dconn.nii
```

## Output Data

Outputs are written to `output_dir` in the same BIDS-like structure:

```
output_dir/
└── sub-<label>/
    └── ses-<label>/
        └── func/
            ├── sub-<label>_ses-<label>_task-<label>_stat-zscored.dconn.nii
            ├── sub-<label>_ses-<label>_task-<label>_ReproTM_template-<template>.dscalar.nii
            ├── sub-<label>_ses-<label>_task-<label>_ReproTM_template-<template>.mat
            ├── sub-<label>_ses-<label>_task-<label>_ReproTM_template-<template>_refine-SCAN.dscalar.nii  # if --refineSCAN
            ├── sub-<label>_ses-<label>_task-<label>_ReproTM_template-<template>_recolored_minsize<N>.dscalar.nii
            └── sub-<label>_ses-<label>_task-<label>_ReproTM_template-<template>_recolored_minsize<N>.dlabel.nii
```


## Usage

```bash
python run.py <bids_dir> <output_dir> participant [options]
```

### Minimal example

```bash
python run.py \
  /data/bids \
  /data/derivatives/repro_tm \
  participant \
  --wb_command /path/to/wb_command \
  --template_infile /path/to/tpl-MyTemplate_seedmaps.mat \
  --template_networks "DMN Vis FP NaN DAN NaN VAN Sal AMN SMd SMl Aud Tpole MTL PMN PON NaN SCAN"
```

### Process specific participants and sessions

```bash
python run.py /data/bids /data/derivatives/repro_tm participant \
  --participant_label sub01 sub02 \
  --session_label ses01 \
  --task_label rest \
  --wb_command /path/to/wb_command \
  --template_infile /path/to/template.mat \
  --template_networks "DMN Vis FP NaN DAN NaN VAN Sal AMN SMd SMl Aud Tpole MTL PMN PON NaN SCAN"
```

### Process all subjects and sessions

If `--participant_label` or `--session_label` are omitted, the pipeline will process all `sub-` and `ses-` directories found in the BIDS input directory. The script logs this behavior explicitly when no filters are provided.

### Dry run (print commands without executing)

```bash
python run.py /data/bids /data/derivatives/repro_tm participant \
  --wb_command /path/to/wb_command \
  --template_infile /path/to/template.mat \
  --template_networks "DMN Vis FP NaN DAN NaN VAN Sal AMN SMd SMl Aud Tpole MTL PMN PON NaN SCAN" \
  --dry_run
```

## Arguments

### Required

| Argument | Description |
|---|---|
| `bids_dir` | Root BIDS input directory |
| `output_dir` | Output/derivatives directory |
| `analysis_level` | Must be `participant` |
| `--wb_command` | Path to the `wb_command` executable |
| `--template_infile` | Path to the template `.mat` file |
| `--template_networks` | Space-separated network names matching template columns |

### Participant Filtering

| Argument | Description |
|---|---|
| `--participant_label` | One or more subject labels (without `sub-` prefix). If omitted, all `sub-` directories in the BIDS root are processed. |
| `--session_label` | One or more session labels (without `ses-` prefix). If omitted, all `ses-` directories under each subject are processed. |
| `--task_label` | One or more task labels to filter dconn files |

### Pipeline Step Toggles

| Argument | Description |
|---|---|
| `--skip_zscore` | Skip step 1: z-scoring |
| `--skip_repro_tm` | Skip step 2: template matching |
| `--skip_minsize` | Skip step 3: min-size cleanup |
| `--skip_dlabel` | Skip step 4: dscalar → dlabel conversion |

When a step is skipped, the pipeline expects its outputs to already exist at the standard output paths.

### Template Matching Options

| Argument | Default | Description |
|---|---|---|
| `--surface_only` | off | Run template matching on cortex only (no subcortex) |
| `--template_thresholding` | off | Threshold template seedmaps to top connections |
| `--template_minthreshold` | `1` | Minimum threshold for template thresholding |
| `--refineSCAN` | off | Re-run matching with higher threshold for SCAN, SMl, SMd |
| `--refineSCAN_minthreshold` | `3` | Threshold applied during SCAN/SMl/SMd refinement |

### Min-size Cleanup Options

| Argument | Default | Description |
|---|---|---|
| `--minsize` | `30` | Minimum network size in greyordinates; network clusters below minimum size are identified and assigned to the mode network assignment of neighboring greyoridinates |

### dscalar → dlabel Conversion Options

| Argument | Default | Description |
|---|---|---|
| `--label_file` | `./dscalar2dlabel/support_files/tpl-ABCC2026-a3-9to16_15networks_label_list.txt` | Label list `.txt` file for `wb_command -cifti-label-import` |

### Miscellaneous

| Argument | Description |
|---|---|
| `--dry_run` | Print commands without executing them |
| `--verbose` | Enable debug-level logging |


## Template Name Inference

The template name embedded in output filenames is extracted automatically from `--template_infile`:

- Files containing `tpl-` → uses the token after `tpl-` up to the next `_`
- Files containing `seedmaps-` → uses the token after `seedmaps-`
- Otherwise → uses the full file stem

## Notes

- The pipeline currently supports only the `participant` analysis level.
- Step 3 (min-size cleanup) uses MATLAB internally. The Python `minsize_v1.0.0.py` script path is also resolved but the MATLAB path (`/common/software/install/migrated/matlab/R2019a/bin/matlab`) is hardcoded; update `step3_minsize_matlab` in `run.py` if your MATLAB installation is elsewhere.
- If `--refineSCAN` is enabled, both the base dscalar and the SCAN-refined dscalar are passed through steps 3 and 4 independently.
