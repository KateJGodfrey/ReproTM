"""
Refactored Python conversion of makeCiftiTemplates_RH.m
Creates CIFTI templates from timeseries data with motion censoring

Key improvements:
- Separated into logical classes and modules
- Improved error handling and validation
- Better memory efficiency
- Vectorized operations where possible
- Type hints for better code clarity
- Reduced code duplication
"""

import numpy as np
import nibabel as nib
from nibabel import cifti2
import os
from pathlib import Path
from scipy import stats
from scipy.io import loadmat, savemat
import random
import time
from typing import List, Tuple, Dict, Optional, Union
from dataclasses import dataclass
import logging
import subprocess

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class Config:
    """Configuration parameters for CIFTI template creation."""
    FD_threshold: float = 0.2
    FD_column: int = 21
    check_motion_first: int = 1
    minutes_to_use: float = 10.0
    ncortgrey: int = 59412
    
    # Region sizes for z-scoring
    L_size: int = 29696
    R_size: int = 29716
    S_size: int = 31870
    
    # Network names
    network_names: List[str] = None
    empty_network_indices: List[int] = None
    
    def __post_init__(self):
        if self.network_names is None:
            self.network_names = ['DMN', 'Vis', 'FP', '', 'DAN', '', 'VAN', 'Sal',
                                 'CO', 'SMd', 'SMl', 'Aud', 'Tpole', 'MTL', 'PMN', 'PON']
        if self.empty_network_indices is None:
            self.empty_network_indices = [3, 5, 16]
    
    def add_scan_network(self):
        """Add SCAN network to network names."""
        if len(self.network_names) < 18:
            self.network_names.extend(['', 'SCAN'])


class PathManager:
    """Manages file paths and settings."""
    
    def __init__(self, support_folder: str = None):
        self.support_folder = support_folder or \
            '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files'
        
        self.paths = {
            'gifti': '/projects/standard/faird/shared/code/external/utilities/gifti-1.6',
            'matlab_cifti': '/projects/standard/faird/shared/code/internal/utilities/Matlab_CIFTI',
            'effect_size': '/projects/standard/faird/shared/code/external/utilities/Matlab_effect_size_toolbox/',
            'consensus_template': f'{self.support_folder}/Networks_template_cleaned.dscalar.nii',
            'consensus_scan': f'{self.support_folder}/Networks_template_cleaned_wABCDscan.dscalar.nii',
            'surface_only_template': f'{self.support_folder}/91282_Greyordinates_surf_only.dtseries.nii',
            'wb_command': '/projects/standard/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'
        }
    
    def get_consensus_path(self, include_scan: bool) -> str:
        """Get path to consensus template based on whether SCAN is included."""
        return self.paths['consensus_scan'] if include_scan else self.paths['consensus_template']


class FileHandler:
    """Handles file I/O operations."""
    
    @staticmethod
    def load_conc_file(filepath: str) -> List[str]:
        """Load concatenated file list."""
        ext = Path(filepath).suffix
        
        if ext == '.conc':
            with open(filepath, 'r') as f:
                return [line.strip() for line in f if line.strip()]
        elif ext == '.mat':
            mat_data = loadmat(filepath)
            key = 'allmasks_outliers_removed_FD02'
            if key in mat_data:
                return [mat_data[key][i, 0] for i in range(mat_data[key].shape[0])]
        
        return [filepath]
    
    @staticmethod
    def validate_files(filepaths: List[str], file_type: str = "file") -> Tuple[List[str], List[int]]:
        """Validate that files exist and return valid files and their indices."""
        valid_files = []
        valid_indices = []
        
        for i, filepath in enumerate(filepaths):
            if isinstance(filepath, (np.ndarray, list)) or Path(filepath).exists():
                valid_files.append(filepath)
                valid_indices.append(i)
            else:
                logger.warning(f"{file_type} {i} does not exist: {filepath}")
        
        logger.info(f"{len(valid_files)} of {len(filepaths)} {file_type}s found")
        return valid_files, valid_indices
    
    @staticmethod
    def load_cifti(filepath: str) -> Dict:
        """Load CIFTI file using nibabel."""
        img = nib.load(filepath)
        return {
            'data': img.get_fdata().squeeze(),
            'header': img.header,
            'img': img
        }
    
    @staticmethod
    def save_cifti(data: np.ndarray, template_cifti: nib.Cifti2Image, 
                   output_path: str, label: str = "Network"):
        """Save CIFTI file with proper structure."""
        # print('Using template cifti for header info:', template_cifti)
        # print('avg data shape:', data.shape)
        # Ensure correct shape (1, n_grayordinates)
        data_reshaped = data.reshape(1, -1)
        # print('reshaped data shape:', data_reshaped.shape)
        
        # Get brain models axis from template
        brain_models_axis = template_cifti.header.get_axis(1)
        # print('template header axes:', brain_models_axis)
        
        # Create scalar axis
        scalar_axis = nib.cifti2.ScalarAxis([label])
        
        # Create new header and image
        new_header = nib.cifti2.Cifti2Header.from_axes((scalar_axis, brain_models_axis))
        new_cifti = nib.Cifti2Image(data_reshaped, new_header)

        # print(f"New CIFTI shape: {new_cifti.shape}")
        # print(f"New CIFTI : {new_cifti}")
        
        nib.save(new_cifti, output_path)
        logger.info(f"Saved CIFTI file to: {output_path}")


class MotionHandler:
    """Handles motion censoring and quality checks."""
    
    def __init__(self, config: Config):
        self.config = config
    
    def load_motion_mask(self, motion_file: Union[str, np.ndarray], 
                        combined_mask: bool = False) -> np.ndarray:
        """Load motion mask from file or array."""
        # Handle array input
        if isinstance(motion_file, (np.ndarray, list)):
            return np.array(motion_file).astype(bool)
        
        # Handle file input
        ext = Path(motion_file).suffix
        
        if ext == '.txt':
            return np.loadtxt(motion_file).astype(bool)
        elif ext == '.mat':
            return self._load_mat_motion(motion_file, combined_mask)
        
        raise ValueError(f"Unsupported motion file format: {ext}")
    
    def _load_mat_motion(self, mat_file: str, combined_mask: bool) -> np.ndarray:
        """Load motion data from .mat file."""
        motion_mat = loadmat(mat_file)
        motion_data = motion_mat['motion_data']
        
        # Find correct FD threshold column
        fd_idx = self._find_fd_column(motion_data)
        
        # Get frame removal vector
        mask_key = 'combined_removal' if combined_mask else 'frame_removal'
        fd_vec = motion_data[0, fd_idx][mask_key]
        # print('fd_vec:',fd_vec)
        
        # Extract the actual array from nested structure
        fd_vec = fd_vec[0, 0]  # Get the inner array
        # print('extracted fd_vec shape:', fd_vec.shape)
        
        # Clean and convert
        fd_vec = self._clean_array(fd_vec).astype(float).flatten()  # Also flatten to ensure 1D
        
        # Convert "1=remove" to "1=keep"
        fd_vec = 1 - fd_vec
        fd_vec = np.clip(fd_vec, 0, 1)
        
        logger.info(f"Motion mask: {np.sum(fd_vec)} good frames of {len(fd_vec)}")
        return fd_vec.astype(bool)
    
    def _find_fd_column(self, motion_data: np.ndarray) -> int:
        """Find column index matching FD threshold."""
        n_cols = motion_data.shape[1]
        all_fd = np.zeros(n_cols)
        
        for i in range(n_cols):
            fd_val = self._extract_scalar(motion_data[0, i]['FD_threshold'])
            all_fd[i] = float(fd_val)
        
        matches = np.where(np.round(all_fd, 3) == np.round(self.config.FD_threshold, 3))[0]
        
        if len(matches) == 0:
            raise ValueError(f"FD threshold {self.config.FD_threshold} not found. "
                           f"Available: {all_fd}")
        
        return matches[0]
    
    @staticmethod
    def _clean_array(arr: np.ndarray) -> np.ndarray:
        """Clean nested arrays and convert to 1D."""
        arr = np.asarray(arr).flatten()
        
        if arr.dtype == object:
            if arr.size == 1:
                return np.array([arr.item()])
            arr = np.concatenate([np.asarray(v).reshape(-1) for v in arr])
        
        return arr
    
    @staticmethod
    def _extract_scalar(value) -> float:
        """Extract scalar value from potentially nested array."""
        while isinstance(value, np.ndarray):
            value = value.item() if value.size == 1 else value.flatten()[0]
        return float(value)
    
    def check_motion_quality(self, motion_files: List[str], TR: float, 
                            project_dir: str, output_name: str) -> Tuple[List[int], List[int]]:
        """Check which subjects meet motion criteria."""
        good_idx = []
        bad_idx = []
        
        for i, motion_file in enumerate(motion_files):
            minutes, _ = self._compute_good_minutes(motion_file, TR)
            logger.info(f"Subject {i}: {minutes:.2f} minutes of good data")
            
            if minutes >= self.config.minutes_to_use:
                good_idx.append(i)
            else:
                bad_idx.append(i)
        
        logger.info(f"{len(good_idx)} subjects pass motion criteria "
                   f"({self.config.minutes_to_use} min at FD={self.config.FD_threshold})")
        
        return good_idx, bad_idx
    
    def _compute_good_minutes(self, motion_file: str, TR: float) -> Tuple[float, float]:
        """Compute minutes of good data after censoring."""
        mask = self.load_motion_mask(motion_file, combined_mask=False)
        good_frames = np.sum(mask)
        good_minutes = (good_frames * TR) / 60
        return good_minutes, 0.0  # mean_FD placeholder


class SubjectMatcher:
    """Matches subjects with their motion files."""
    
    @staticmethod
    def match_subjects(timeseries_files: List[str], 
                       motion_files: List[str]) -> Tuple[List[str], List[str]]:
        """Match timeseries files with corresponding motion files."""
        matched_ts = []
        matched_motion = []
        
        for ts_file in timeseries_files:
            ts_id = SubjectMatcher._extract_subject_id(ts_file)
            
            for motion_file in motion_files:
                if isinstance(motion_file, (np.ndarray, list)):
                    continue
                
                motion_id = SubjectMatcher._extract_subject_id(motion_file)
                
                if ts_id == motion_id:
                    matched_ts.append(ts_file)
                    matched_motion.append(motion_file)
                    logger.debug(f"Matched subject: {ts_id}")
                    break
        
        logger.info(f"Matched {len(matched_ts)} subjects with motion files")
        return matched_ts, matched_motion
    
    @staticmethod
    def _extract_subject_id(filepath: str) -> str:
        """Extract subject ID from filepath."""
        parts = filepath.split('/')
        try:
            func_idx = parts.index('func')
            return '/'.join(parts[:func_idx]) if func_idx > 0 else parts[0]
        except ValueError:
            return parts[0]


class NetworkAnalyzer:
    """Computes network connectivity maps."""
    
    def __init__(self, config: Config, consensus_data: np.ndarray):
        self.config = config
        self.consensus_data = consensus_data
        # print('consensus_data shape:', consensus_data.shape)
    
    def compute_subject_seedmaps(self, timeseries_data: np.ndarray, 
                                 surface_only: bool = False) -> np.ndarray:
        """Compute correlation maps for all networks for a single subject."""
        n_grayordinates = timeseries_data.shape[0]
        n_networks = len(self.config.network_names)
        seedmaps = np.zeros((n_grayordinates, n_networks))
        
        for net_idx, net_name in enumerate(self.config.network_names):
            if net_idx in self.config.empty_network_indices:
                continue
            
            seedmaps[:, net_idx] = self._compute_network_correlation(
                timeseries_data, net_idx + 1, surface_only
            )
        
        return seedmaps
    
    def _compute_network_correlation(self, timeseries_data: np.ndarray, 
                                    net_label: int, surface_only: bool) -> np.ndarray:
        """Compute correlation map for a single network."""
        # Get network voxels
        net_mask = self.consensus_data == net_label
        net_idx = np.where(net_mask)[0]
        
        if len(net_idx) == 0:
            logger.warning(f"No voxels found for network {net_label}")
            return np.zeros(timeseries_data.shape[0])
        
        # Compute network average timeseries
        net_avg = np.nanmean(timeseries_data[net_idx, :], axis=0)
        
        # Compute correlations
        n_voxels = self.config.ncortgrey if surface_only else timeseries_data.shape[0]
        corr_map = np.zeros(timeseries_data.shape[0])
        
        for vox in range(n_voxels):
            corr_map[vox] = self._pearson_correlation(net_avg, timeseries_data[vox, :])
        
        return corr_map
    
    @staticmethod
    def _pearson_correlation(x: np.ndarray, y: np.ndarray) -> float:
        """Compute Pearson correlation, handling NaNs."""
        valid = ~(np.isnan(x) | np.isnan(y))
        
        if np.sum(valid) < 2:
            return np.nan
        
        return np.corrcoef(x[valid], y[valid])[0, 1]
    
    def aggregate_seedmaps(self, subject_seedmaps: List[np.ndarray], 
                          zscore_regions: bool = False,
                          surface_only: bool = False) -> np.ndarray:
        """Aggregate seedmaps across subjects using Fisher z-transform."""
        n_grayordinates = subject_seedmaps[0].shape[0]
        n_networks = len(self.config.network_names)
        avg_seedmaps = np.zeros((n_grayordinates, n_networks))
        
        for net_idx in range(n_networks):
            if net_idx in self.config.empty_network_indices:
                continue
            
            # Fisher z-transform, average, inverse transform
            z_maps = [np.arctanh(sm[:, net_idx]) for sm in subject_seedmaps]
            z_avg = np.mean(z_maps, axis=0)
            avg_seedmaps[:, net_idx] = np.tanh(z_avg)
            
            # Z-score within regions if requested
            if zscore_regions:
                avg_seedmaps[:, net_idx] = self._zscore_regions(
                    avg_seedmaps[:, net_idx], surface_only
                )
        
        return avg_seedmaps
    
    def _zscore_regions(self, data: np.ndarray, surface_only: bool) -> np.ndarray:
        """Z-score within left, right, and subcortical regions."""
        result = data.copy()
        
        # Left hemisphere
        result[:self.config.L_size] = stats.zscore(result[:self.config.L_size])
        
        # Right hemisphere
        result[self.config.L_size:self.config.ncortgrey] = \
            stats.zscore(result[self.config.L_size:self.config.ncortgrey])
        
        # Subcortical (if not surface only)
        if not surface_only and len(result) > self.config.ncortgrey:
            result[self.config.ncortgrey:] = stats.zscore(result[self.config.ncortgrey:])
        
        return result


class CIFTITemplateCreator:
    """Main class for creating CIFTI templates."""
    
    def __init__(self, config: Config, path_manager: PathManager):
        self.config = config
        self.paths = path_manager
        self.file_handler = FileHandler()
        self.motion_handler = MotionHandler(config)
        self.subject_matcher = SubjectMatcher()
    
    def create_templates(self, timeseries_conc: str, motion_conc: str,
                        project_dir: str, TR: float,
                        zscore_regions: bool = False,
                        power_motion: bool = True,
                        remove_outlier: bool = True,
                        surface_only: bool = False,
                        use_motion_criteria: bool = True,
                        combined_mask: bool = False,
                        include_scan: bool = False):
        """Main pipeline for creating CIFTI templates."""
        
        # Setup
        os.makedirs(project_dir, exist_ok=True)
        
        if include_scan:
            self.config.add_scan_network()
        
        # Load file lists
        ts_files = self.file_handler.load_conc_file(timeseries_conc)
        motion_files = self.file_handler.load_conc_file(motion_conc)
        
        # Validate files
        ts_files, _ = self.file_handler.validate_files(ts_files, "timeseries")
        motion_files, _ = self.file_handler.validate_files(motion_files, "motion")
        
        # Match subjects
        ts_files, motion_files = self.subject_matcher.match_subjects(ts_files, motion_files)
        
        if len(ts_files) != len(motion_files):
            raise ValueError(f"Mismatch: {len(ts_files)} timeseries, "
                           f"{len(motion_files)} motion files")
        
        # Check motion quality
        if self.config.check_motion_first:
            good_idx, bad_idx = self.motion_handler.check_motion_quality(
                motion_files, TR, project_dir, Path(motion_conc).stem
            )
            
            if use_motion_criteria:
                ts_files = [ts_files[i] for i in good_idx]
                motion_files = [motion_files[i] for i in good_idx]
        
        # Load consensus template
        consensus_path = self.paths.get_consensus_path(include_scan)
        print('paths:',self.paths.paths['wb_command'])
        # print('consensus_path:',consensus_path)
        consensus = self.file_handler.load_cifti(consensus_path)
        
        # Create analyzer
        analyzer = NetworkAnalyzer(self.config, consensus['data'])
        
        # Process subjects or load cached
        seedmaps_file = Path(project_dir) / f"seedmaps_{Path(timeseries_conc).stem}.mat"
        
        if seedmaps_file.exists():
            logger.info("Loading cached seedmaps")
            subject_seedmaps = self._load_cached_seedmaps(seedmaps_file)
        else:
            logger.info(f"Processing {len(ts_files)} subjects")
            subject_seedmaps, masks_after, masks_before = self._process_all_subjects(
                ts_files, motion_files, analyzer, surface_only, 
                power_motion, combined_mask, TR, project_dir, timeseries_conc, remove_outlier
            )
            self._save_seedmaps(seedmaps_file, subject_seedmaps, masks_after, masks_before)
        
        # Remove subjects with NaNs
        clean_seedmaps, clean_files = self._remove_nan_subjects(
            subject_seedmaps, ts_files
        )
        
        # Aggregate and save results
        avg_seedmaps = analyzer.aggregate_seedmaps(
            clean_seedmaps, zscore_regions, surface_only
        )
        
        # self._save_results(avg_seedmaps, clean_files, project_dir, 
        #                   Path(timeseries_conc).stem, consensus['img'],
        #                   zscore_regions, surface_only)

        ######## Using first subject as template to save CIFTI files ########
        self._save_results(avg_seedmaps, clean_files, motion_files, project_dir, 
                          Path(timeseries_conc).stem,
                          zscore_regions, surface_only, combined_mask, 
                          self.config.FD_threshold, include_scan, self.config.minutes_to_use, 
                          power_motion, remove_outlier,ts_files, use_motion_criteria, good_idx, bad_idx)
        
        logger.info("Template creation complete!")
    
    def _process_all_subjects(self, ts_files: List[str], motion_files: List[str],
                             analyzer: NetworkAnalyzer, surface_only: bool,
                             power_motion: bool, combined_mask: bool,
                             TR: float, project_dir: str, timeseries_conc: str, remove_outliers: bool) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
        """Process all subjects to compute seedmaps."""
        subject_seedmaps = []
        allmasks_outliers_removed_FD02 = []
        allmasks_before_outliers_removed_FD02 = []

        for i, (ts_file, motion_file) in enumerate(zip(ts_files, motion_files)):
            start = time.time()
            logger.info(f"Processing subject {i+1}/{len(ts_files)}: {ts_file}")
            orig_motion_filename = os.path.splitext(os.path.basename(motion_file))[0]
            
            # Load timeseries
            ts = self.file_handler.load_cifti(ts_file)
            ts_data = ts['data']
            
            # Ensure correct orientation (grayordinates x timepoints)
            if ts_data.shape[0] < ts_data.shape[1]:
                ts_data = ts_data.T
            
            # Load and apply motion mask
            motion_mask = self.motion_handler.load_motion_mask(
                motion_file, combined_mask
            )
            allmasks_before_outliers_removed_FD02.append(motion_mask)

            # ========== Outlier detection ==========
            if combined_mask == False:
                if remove_outliers == True:
                    logger.info("Removal outliers not specified. It will be performed by default.")
                    stdev_temp_filename = f"{Path(timeseries_conc).stem}_temp.txt"
                    motion_mask = self.CensorBOLDoutliers(self.paths.paths['wb_command'], ts_files, i, stdev_temp_filename, motion_mask)
                else:
                    logger.info("Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed.")
                
                good_frames_idx = np.where(motion_mask == 1)[0]
                good_minutes = (len(good_frames_idx) * TR) / 60
                
                if good_minutes < 0.5:
                    logger.warning(f"Subject {i} has less than 30 seconds of good data")
                    continue
                elif self.config.minutes_to_use > good_minutes:
                    # Use all available frames
                    output_file = os.path.join(project_dir, 
                        f"{orig_motion_filename}_{self.config.FD_threshold}_cifti_censor_FD_vector_All_Good_Frames.txt")
                    np.savetxt(output_file, motion_mask, fmt='%1.0f')
                else:
                    # Randomly sample frames to match minutes_to_use
                    good_frames_needed = round(self.config.minutes_to_use * 60 / TR)
                    rand_good_frames = sorted(random.sample(list(range(len(good_frames_idx))), good_frames_needed))
                    motion_mask_cut = np.zeros(len(motion_mask))
                    ones_idx = good_frames_idx[rand_good_frames]
                    motion_mask_cut[ones_idx] = 1
                    
                    output_file = os.path.join(project_dir, 
                        f"{orig_motion_filename}_{self.config.FD_threshold}_cifti_censor_FD_vector_{self.config.minutes_to_use}_minutes_of_data_at_{self.config.FD_threshold}_threshold.txt")
                    np.savetxt(output_file, motion_mask_cut, fmt='%1.0f')
                    motion_mask = motion_mask_cut
            else:
                logger.info("Frames with outliers in BOLD std dev not removed. Maybe you have already performed outlier detection.")
            
            allmasks_outliers_removed_FD02.append(motion_mask)

            # Handle mask/timeseries length mismatch
            if len(motion_mask) != ts_data.shape[1]:
                logger.warning(f"Mask length ({len(motion_mask)}) != "
                             f"timepoints ({ts_data.shape[1]})")
                if len(motion_mask) < ts_data.shape[1]:
                    logger.error("Mask too short, skipping subject")
                    continue
                motion_mask = motion_mask[:ts_data.shape[1]]

            # Censor timeseries
            ts_data = ts_data[:, motion_mask]
            logger.info(f"Censored timeseries: {ts_data.shape}")
            
            # Compute seedmaps
            seedmaps = analyzer.compute_subject_seedmaps(ts_data, surface_only)
            subject_seedmaps.append(seedmaps)
            
            logger.info(f"Subject {i+1} completed in {time.time()-start:.2f}s")
        
        return subject_seedmaps, allmasks_outliers_removed_FD02, allmasks_before_outliers_removed_FD02
    
    def isthisanoutlier(data, method='median'):
        """
        Detect outliers using the median absolute deviation method.
        
        Parameters:
        -----------
        data : array-like
            Input data array
        method : str
            Method for outlier detection ('median')
            
        Returns:
        --------
        ndarray : Boolean array where True indicates an outlier
        """
        if method == 'median':
            # Calculate median absolute deviation
            median = np.median(data)
            mad = np.median(np.abs(data - median))
            
            # Modified Z-score method (similar to MATLAB's isoutlier with 'median')
            # Using threshold of 3 (common default)
            if mad == 0:
                # If MAD is 0, use standard deviation
                std = np.std(data)
                if std == 0:
                    return np.zeros(len(data), dtype=bool)
                modified_z_scores = 0.6745 * (data - median) / std
            else:
                modified_z_scores = 0.6745 * (data - median) / mad
            
            return np.abs(modified_z_scores) > 3
        else:
            raise ValueError(f"Method '{method}' not implemented")


    def CensorBOLDoutliers(self, wb_command, A, i, stdev_temp_filename, FDvec):
        """
        Provides an additional level of data cleaning to correlations made from 
        the dtseries beyond FD censoring.
        
        This code calls workbench command to calculate the standard deviation
        of the dtseries. Motion censoring appears to remove large artifacts
        introduced into the magnetic field. However, some artifacts are present 
        in the BOLD data even after frames >0.2 (or whichever threshold) are removed.
        These frames can be identified by their unusually high standard deviation.
        
        Parameters:
        -----------
        wb_command : str
            Path to wb_command executable
        A : str or list
            Subject's dtseries file path (string) or list of dtseries files (list)
        i : int
            Which subject in the list you're currently using (only used if A is a list)
        stdev_temp_filename : str
            Name of temporary file that wb_command will write standard deviations to
        FDvec : ndarray
            Motion censored vector (1s and 0s). Should already be motion censored 
            and have the same number of values (frames) as your dtseries.
            
        Returns:
        --------
        FDvec : ndarray
            New vector of 1s and 0s with outliers removed
            
        Author: R. Hermosillo
        Original MATLAB version: 6/20/2018
        Updated: 9/24/2019
        Python conversion: 2026
        """
        
        # Construct the command based on input type
        if isinstance(A, list):
            cmd = f"{wb_command} -cifti-stats {A[i]} -reduce STDEV > {stdev_temp_filename}"
        elif isinstance(A, str):
            cmd = f"{wb_command} -cifti-stats {A} -reduce STDEV > {stdev_temp_filename}"
        else:
            raise TypeError('Input time series name must be a string or list.')
        
        # Execute the command
        subprocess.run(cmd, shell=True, check=True)
        
        # Give time to write temp file
        time.sleep(3)
        
        # Try to load the file with multiple retries
        max_attempts = 4
        wait_times = [3, 10, 10, 10]  # seconds to wait between attempts
        
        STDEV_file = None
        for attempt in range(max_attempts):
            try:
                STDEV_file = np.loadtxt(stdev_temp_filename)
                print('STDEV file loaded successfully.')
                break
            except Exception as e:
                if attempt < max_attempts - 1:
                    attempt_num = attempt + 1
                    print(f'Attempt to load stdev file failed: {stdev_temp_filename}')
                    print(f'System may be slow- retrying in {wait_times[attempt+1]} seconds ({attempt_num}/{max_attempts-1})...')
                    time.sleep(wait_times[attempt + 1])
                else:
                    print(f'Unable to load stdev file: {stdev_temp_filename}')
                    print('Exiting.')
                    return FDvec
        
        # Clean up temporary file
        try:
            os.remove(stdev_temp_filename)
        except OSError as e:
            print(f"Warning: Could not remove temporary file {stdev_temp_filename}: {e}")
        
        # Find the kept frames from the FD mask
        FDvec_keep_idx = np.where(FDvec == 1)[0]
        
        # Find outliers using the median method
        Outlier_file = self.isthisanoutlier(STDEV_file[FDvec_keep_idx], 'median')
        
        # Find outlier indices
        Outlier_idx = np.where(Outlier_file == 1)[0]
        
        # Set outliers to zero within FDvec
        FDvec[FDvec_keep_idx[Outlier_idx]] = 0
        
        return FDvec
    
    def _load_cached_seedmaps(self, filepath: Path) -> List[np.ndarray]:
        """Load cached seedmaps from file."""
        mat_data = loadmat(str(filepath))
        return [mat_data['seedmapstimeseries'][i, 0] 
                for i in range(mat_data['seedmapstimeseries'].shape[0])]
    
    def _save_seedmaps(self, filepath: Path, seedmaps: List[np.ndarray], 
                   masks_after: List[np.ndarray], masks_before: List[np.ndarray]):
        """Save seedmaps and motion masks to file.
        
        Args:
            filepath: Path to save the .mat file
            seedmaps: List of seedmap arrays for each subject
            masks_after: Motion masks after outlier removal
            masks_before: Original motion masks before outlier removal
        """
        # Convert seedmaps to cell array format
        seedmaps_cell = np.empty((len(seedmaps), 1), dtype=object)
        for i, sm in enumerate(seedmaps):
            seedmaps_cell[i, 0] = sm
        
        # Convert masks to cell array format
        masks_after_cell = np.empty((len(masks_after), 1), dtype=object)
        for i, mask in enumerate(masks_after):
            masks_after_cell[i, 0] = mask
        
        masks_before_cell = np.empty((len(masks_before), 1), dtype=object)
        for i, mask in enumerate(masks_before):
            masks_before_cell[i, 0] = mask
        
        # Save all data to mat file
        savemat(str(filepath), {
            'seedmapstimeseries': seedmaps_cell,
            'allmasks_outliers_removed_FD02': masks_after_cell,
            'allmasks_before_outliers_removed_FD02': masks_before_cell
        })
        
        logger.info(f"Saved seedmaps and motion masks to {filepath}")
    
    def _remove_nan_subjects(self, subject_seedmaps: List[np.ndarray],
                            ts_files: List[str]) -> Tuple[List[np.ndarray], List[str]]:
        """Remove subjects with NaN values."""
        clean_seedmaps = []
        clean_files = []
        
        for i, (seedmap, ts_file) in enumerate(zip(subject_seedmaps, ts_files)):
            has_nan = False
            
            for net_idx in range(len(self.config.network_names)):
                if net_idx in self.config.empty_network_indices:
                    continue
                
                if np.any(np.isnan(seedmap[:, net_idx])):
                    logger.warning(f"Subject {i} has NaNs in network {net_idx}")
                    has_nan = True
                    break
            
            if not has_nan:
                clean_seedmaps.append(seedmap)
                clean_files.append(ts_file)
        
        logger.info(f"Removed {len(subject_seedmaps) - len(clean_seedmaps)} "
                   f"subjects with NaNs")
        return clean_seedmaps, clean_files
    
    def _save_results(self, avg_seedmaps: np.ndarray, clean_files: List[str],
                 motion_files: List[str], project_dir: str, base_name: str,
                 zscore_regions: bool, surface_only: bool, combined_mask: bool,
                 FD_threshold: float, include_scan: bool, minutes_to_use: float,
                 power_motion: bool, remove_outliers: bool, ts_files: List[str],
                 use_motion_criteria: bool, good_idx: List[int], bad_idx: List[int]):
        """Save final results including seedmaps and metadata."""
        suffix = "_Zscored" if zscore_regions else ""
        n_subjects = len(clean_files)

        template_cifti = nib.load(clean_files[0])
        
        # Save individual network maps
        for net_idx, net_name in enumerate(self.config.network_names):
            if net_idx in self.config.empty_network_indices or not net_name:
                continue
            
            output_path = Path(project_dir) / \
                f"seedmaps_{base_name}_{net_name}_network{suffix}.dtseries.nii"
            
            self.file_handler.save_cifti(
                avg_seedmaps[:, net_idx], template_cifti, str(output_path), net_name
            )
        
        # Save comprehensive summary .mat file
        summary_path = Path(project_dir) / \
            f"seedmaps_{base_name}_n_{n_subjects}_all_networks{suffix}.mat"
        
        savemat(str(summary_path), {
            'B': np.array(motion_files, dtype=object).reshape(1, -1),
            'bad_subs_idx': np.array(bad_idx, dtype=np.float64).reshape(1, -1),
            'cleansubs': np.array(clean_files, dtype=object).reshape(1, -1),
            'combined_outliermask_provided': combined_mask,
            'FD_threshold': FD_threshold,
            'good_subs_idx': np.array(good_idx, dtype=np.float64).reshape(1, -1),
            'include_scan_net': include_scan,
            'minutes_to_use': minutes_to_use,
            'power_motion': power_motion,
            'remove_outliers': remove_outliers,
            'seed_matrix': avg_seedmaps,
            'subs': np.array(ts_files, dtype=object).reshape(1, -1),
            'surface_only': surface_only,
            'use_only_subjects_that_pass_motion_criteria': use_motion_criteria,
            'Zscore_regions': zscore_regions
            # 'network_names': np.array(self.config.network_names, dtype=object),
        })
        
        logger.info(f"Saved results for {n_subjects} subjects to {project_dir}")


def main():
    """Main entry point with argument parsing."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Create CIFTI templates with motion censoring'
    )
    parser.add_argument('--timeseries', required=True, 
                       help='Path to timeseries conc file')
    parser.add_argument('--TR', type=float, required=True, 
                       help='Repetition time in seconds')
    parser.add_argument('--motion', required=True, 
                       help='Path to motion conc file')
    parser.add_argument('--project_dir', required=True, 
                       help='Output directory')
    parser.add_argument('--zscore_regions', action='store_true',
                       help='Z-score within regions')
    parser.add_argument('--power_motion', action='store_true', default=True,
                       help='Use Power motion method')
    parser.add_argument('--remove_outliers', action='store_true', default=True,
                       help='Remove outliers based on BOLD std dev')
    parser.add_argument('--surface_only', action='store_true',
                       help='Surface only')
    parser.add_argument('--use_motion_criteria', action='store_true', default=True,
                       help='Filter by motion criteria')
    parser.add_argument('--combined_mask', action='store_true',
                       help='Combined outlier mask provided')
    parser.add_argument('--include_scan', action='store_true',
                       help='Include SCAN network')
    
    args = parser.parse_args()
    
    # Initialize components
    config = Config()
    path_manager = PathManager()
    creator = CIFTITemplateCreator(config, path_manager)
    
    # Run pipeline
    creator.create_templates(
        timeseries_conc=args.timeseries,
        motion_conc=args.motion,
        project_dir=args.project_dir,
        TR=args.TR,
        zscore_regions=args.zscore_regions,
        power_motion=args.power_motion,
        remove_outlier=args.remove_outliers,
        surface_only=args.surface_only,
        use_motion_criteria=args.use_motion_criteria,
        combined_mask=args.combined_mask,
        include_scan=args.include_scan
    )


if __name__ == "__main__":
    main()