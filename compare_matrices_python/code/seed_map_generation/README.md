# CIFTI Template Creator

Creates network connectivity templates from fMRI timeseries data with motion censoring.

## Installation

```bash
pip install numpy scipy nibabel
```

## Quick Start

```bash
python makeCiftiTemplate_RH_refactor.py \
  --timeseries /path/to/timeseries.conc \
  --motion /path/to/motion.conc \
  --project_dir /path/to/output \
  --TR 0.8
```

## Required Arguments

| Flag | Description |
|------|-------------|
| `--timeseries` | Path to timeseries file list (.conc, .mat, or single CIFTI file) |
| `--motion` | Path to motion data file list (.conc, .mat, or .txt with binary masks) |
| `--project_dir` | Output directory for results |
| `--TR` | Repetition time in seconds (e.g., 0.8) |

## Optional Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--zscore_regions` | Off | Z-score maps within left/right hemisphere and subcortical regions |
| `--power_motion` | On | Use Power et al. FD-based motion censoring (FD < 0.2mm) |
| `--remove_outliers` | On | Remove frames with BOLD signal outliers (in addition to FD) |
| `--surface_only` | Off | Process only cortical surface (59,412 grayordinates) |
| `--use_motion_criteria` | On | Exclude subjects with < 10 min of clean data |
| `--combined_mask` | Off | Motion files already include FD + outlier masks |
| `--include_scan` | Off | Add SCAN network (18 networks instead of 16) |

## Common Use Cases

**Standard analysis:**
```bash
python cifti_template_creator.py \
  --timeseries data.conc \
  --motion motion.conc \
  --project_dir results \
  --TR 0.8 \
  --zscore_regions
```

**Cortical-only with all quality filters:**
```bash
python cifti_template_creator.py \
  --timeseries data.conc \
  --motion motion.conc \
  --project_dir results \
  --TR 0.8 \
  --surface_only \
  --use_motion_criteria \
  --remove_outliers
```

**Pre-computed masks (skip outlier detection):**
```bash
python cifti_template_creator.py \
  --timeseries data.conc \
  --motion combined_masks.conc \
  --project_dir results \
  --TR 0.8 \
  --combined_mask
```

## Output Files

- `seedmaps_{name}_{Network}_network.dtseries.nii` - Individual network correlation maps
- `seedmaps_{name}_n_{N}_all_networks.mat` - Summary file with all networks and metadata
- `{motion}_{FD}_cifti_censor_FD_vector_*.txt` - Binary censoring vectors per subject
- `seedmaps_{name}.mat` - Cached results (delete to force recomputation)

## Networks

**Default (14):** DMN, Vis, FP, DAN, VAN, Sal, CO, SMd, SMl, Aud, Tpole, MTL, PMN, PON

**With --include_scan (15):** Above + SCAN

## Configuration

Edit `Config` class to change:
- `FD_threshold` (default: 0.2mm)
- `minutes_to_use` (default: 10.0)
- Network definitions

## Workflow

1. Load and validate timeseries/motion files
2. Match subjects between datasets
3. Check motion quality (exclude subjects with insufficient clean data)
4. Detect BOLD outliers (if enabled)
5. Censor frames based on FD + outliers
6. Compute network correlation maps per subject
7. Aggregate using Fisher z-transform
8. Save CIFTI outputs

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "Mismatch: X timeseries, Y motion" | Ensure matching subject lists in both .conc files |
| "Less than 30 seconds of good data" | Check motion file quality or adjust FD threshold |
| "Mask length != timepoints" | Verify motion data matches timeseries scan |
| Need to reprocess | Delete `seedmaps_{name}.mat` cache file |

## Tips

- **Motion threshold:** 0.2mm standard, 0.15mm for high-motion populations
- **Data minimum:** 10 minutes clean data recommended
- **Z-scoring:** Use `--zscore_regions` for cross-region comparisons
- **Surface-only:** Faster processing, good for cortical network analysis
- **Caching:** Keep .mat cache when only changing aggregation parameters