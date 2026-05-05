"""
Python conversion of makeCiftiTemplates_RH.m
Creates CIFTI templates from timeseries data with motion censoring
"""

import numpy as np
import nibabel as nib
from nibabel import cifti2
import os
import sys
from pathlib import Path
import subprocess
import warnings
from scipy import stats
from scipy.io import loadmat, savemat
import random
import time


def makeCiftiTemplates_RH(dt_or_ptseries_conc_file, TR, all_motion_conc_file, 
                          project_dir, Zscore_regions, power_motion, remove_outliers,
                          surface_only, use_only_subjects_that_pass_motion_criteria,
                          combined_outliermask_provided, include_scan_net):
    """
    Create CIFTI templates from timeseries data with motion censoring.
    
    Parameters:
    -----------
    dt_or_ptseries_conc_file : str
        Path to concatenated timeseries file (.conc) or single file
    TR : float
        Repetition time in seconds
    all_motion_conc_file : str
        Path to motion file(s) (.conc or .mat)
    project_dir : str
        Output directory for results
    Zscore_regions : int
        Whether to z-score within regions (0 or 1)
    power_motion : int
        Use Power method for motion censoring (0 or 1)
    remove_outliers : int
        Remove outlier frames (0 or 1)
    surface_only : int
        Use surface data only (0 or 1)
    use_only_subjects_that_pass_motion_criteria : int
        Filter subjects by motion (0 or 1)
    combined_outliermask_provided : int
        Whether combined outlier mask is provided (0 or 1)
    include_scan_net : int
        Include SCAN network (0 or 1)
    """
    
    # ========== Hardcoded parameters ==========
    FD_threshold = 0.2
    FD_column = 21
    check_motion_first = 1
    minutes_to_use = 10
    
    if surface_only == 1:
        Zscore_regions = 0
    
    if Zscore_regions == 1:
        L_size = 29696  # number of parcellations
        R_size = 29716
        S_size = 31870
    
    ncortgrey = 59412
    
    # ========== Setup paths ==========
    # this_code = os.path.abspath(__file__)
    # code_dir = os.path.dirname(this_code)
    support_folder = '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files'
    
    # Add paths (Python equivalent of MATLAB addpath)
    sys.path.insert(0, support_folder)
    sys.path.insert(0, '/projects/standard/faird/shared/code/external/utilities/MSCcodebase-master/Utilities/')
    
    # Load settings
    settings = settings_comparematrices()
    print('settings:')
    print(settings)
    
    # wb_command path
    wb_command = settings.get('path_wb_c', 'wb_command')
    
    # Network names
    if include_scan_net == 1:
        network_names = ['DMN', 'Vis', 'FP', '', 'DAN', '', 'VAN', 'Sal', 
                        'CO', 'SMd', 'SMl', 'Aud', 'Tpole', 'MTL', 'PMN', 'PON', '', 'SCAN']
    else:
        network_names = ['DMN', 'Vis', 'FP', '', 'DAN', '', 'VAN', 'Sal',
                        'CO', 'SMd', 'SMl', 'Aud', 'Tpole', 'MTL', 'PMN', 'PON']
    
    # ========== Check timeseries files ==========
    conc_ext = os.path.splitext(dt_or_ptseries_conc_file)[1]
    if conc_ext == '.conc':
        with open(dt_or_ptseries_conc_file, 'r') as f:
            subs = [line.strip() for line in f.readlines() if line.strip()]
    else:
        subs = [dt_or_ptseries_conc_file]
    
    # Verify timeseries files exist
    subsfound = 0
    for i, sub in enumerate(subs):
        if not os.path.exists(sub):
            print(f"Subject Series {i} does not exist")
            print(sub)
            continue
        else:
            subsfound += 1
    
    print(f"{subsfound} of {len(subs)} timeseries files found. All series files exist continuing...")
    
    # ========== Load motion files ==========
    motion_ext = os.path.splitext(all_motion_conc_file)[1]
    if motion_ext == '.conc':
        with open(all_motion_conc_file, 'r') as f:
            B = [line.strip() for line in f.readlines() if line.strip()]
        # print(B)
    elif motion_ext == '.mat':
        mat_data = loadmat(all_motion_conc_file)
        B = [mat_data['allmasks_outliers_removed_FD02'][i, 0] for i in range(mat_data['allmasks_outliers_removed_FD02'].shape[0])]
    else:
        B = [all_motion_conc_file]
    
    # ========== Verify motion files exist ==========
    subsfound_motion = 0
    good_subs_list = []
    for i, motion_file in enumerate(B):
        if isinstance(motion_file, (np.ndarray, list)):
            subsfound_motion += 1
        else:
            if not os.path.exists(motion_file):
                print(f"Motion file {i} does not exist: {motion_file}")
                good_subs_list.append(0)
            else:
                subsfound_motion += 1
    
    print(f"{subsfound_motion} of {len(B)} motion files found. Motion files exist continuing...")
    
    # ========== Match subjects with motion files ==========
    good_sub = []
    good_B = []
    
    for i in range(len(subs)):
        if not os.path.exists(subs[i]):
            print(f"Subject Series {i} does not exist")
            print(subs[i])
            continue
        
        # Extract subject ID from timeseries file path
        subs_parts = subs[i].split('/')
        try:
            subs_idx = subs_parts.index('func')
            if subs_idx > 0:
                subs_id = '/'.join(subs_parts[:subs_idx])
            else:
                subs_id = subs_parts[0]
        except ValueError:
            subs_id = subs_parts[0]
        
        # Search for matching subject in motion files
        for j in range(len(B)):
            if isinstance(B[j], (np.ndarray, list)):
                B_parts = []
                B_id = ''
            else:
                B_parts = B[j].split('/')
                try:
                    B_idx = B_parts.index('func')
                    if B_idx > 0:
                        B_id = '/'.join(B_parts[:B_idx])
                    else:
                        B_id = B_parts[0]
                except ValueError:
                    B_id = B_parts[0] if B_parts else ''
            
            # Check if IDs match and motion file exists
            motion_exists = isinstance(B[j], (np.ndarray, list)) or os.path.exists(B[j])
            if subs_id == B_id and motion_exists:
                print('IDs match and motion file exists, adding to good_sub and good_B')
                good_sub.append(subs[i])
                good_B.append(B[j])
                break
    
    print(f"{len(good_sub)} of {len(subs)} subjects matched with motion files.")
    print(f"{len(good_B)} of {len(B)} motion files matched.")
    
    subs = good_sub
    B = good_B
    
    # ========== Verify matching lengths ==========
    if len(subs) == len(B):
        print("Length of motion conc file matches length of series conc file. Good job.")
    else:
        raise ValueError(f"Length of motion conc file ({len(B)}) does not match length of series conc file ({len(subs)}). Check your inputs.")
    
    # ========== Parse file names ==========
    file_name = os.path.basename(dt_or_ptseries_conc_file)
    file_root_no_ext = os.path.splitext(file_name)[0]
    
    if surface_only == 1:
        file_root_no_ext = f"{file_root_no_ext}_SurfOnly"
    
    # ========== Load consensus network template ==========
    if include_scan_net == 1:
        consen_path = '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned_wABCDscan.dscalar.nii'
    else:
        consen_path = settings['path'][5]  # Adjust index as needed
    
    consen = ft_read_cifti_mod(consen_path)
    consen_data = consen['data']
    # print(f"Consensus template shape: {consen_data}")
    
    # ========== Check motion quality ==========
    if check_motion_first == 1:
        output_name = os.path.splitext(os.path.basename(all_motion_conc_file))[0]
        
        good_subs_idx = []
        bad_subs_idx = []
        really_bad_subs_idx = []
        
        for i in range(len(B)):
            all_subjects_minutes, all_mean_FD = check_twins_motion(
                B[i], TR, project_dir, output_name, 
                FD_column, None, 0, 1, 0, 0
            )
            print(f"Subject {i} has {all_subjects_minutes} minutes of data after motion censoring at FD {FD_threshold}")
            
            if minutes_to_use <= all_subjects_minutes:
                good_subs_idx.append(i)
            else:
                bad_subs_idx.append(i)
            
            if all_subjects_minutes < 0.5:
                really_bad_subs_idx.append(i)
        
        print(f"{len(good_subs_idx)} subjects have at least {minutes_to_use} minutes of data after motion censoring at FD {FD_threshold}")
        print(f"{len(bad_subs_idx)} subjects have less than {minutes_to_use} minutes of data after motion censoring at FD {FD_threshold}")
        print(f"{len(really_bad_subs_idx)} subjects have less than 30 seconds of data after motion censoring at FD {FD_threshold}")
    else:
        good_subs_idx = list(range(len(subs)))
    
    # ========== Filter subjects by motion criteria ==========
    if use_only_subjects_that_pass_motion_criteria == 1:
        subs = [subs[i] for i in good_subs_idx]
        B = [B[i] for i in good_subs_idx]
    
    # ========== Load or compute seedmaps ==========
    seedmaps_file = os.path.join(project_dir, f'seedmaps_{file_root_no_ext}.mat')
    if os.path.exists(seedmaps_file):
        print("Loading previously generated data from each subject for template")
        mat_data = loadmat(seedmaps_file)
        seedmapstimeseries = [mat_data['seedmapstimeseries'][i, 0] for i in range(mat_data['seedmapstimeseries'].shape[0])]
        allmasks_outliers_removed_FD02 = [mat_data['allmasks_outliers_removed_FD02'][i, 0] for i in range(mat_data['allmasks_outliers_removed_FD02'].shape[0])]
        allmasks_before_outliers_removed_FD02 = []
        if 'allmasks_before_outliers_removed_FD02' in mat_data:
            allmasks_before_outliers_removed_FD02 = [mat_data['allmasks_before_outliers_removed_FD02'][i, 0] for i in range(mat_data['allmasks_before_outliers_removed_FD02'].shape[0])]
    else:
        seedmapstimeseries = []
        allmasks_outliers_removed_FD02 = []
        allmasks_before_outliers_removed_FD02 = []
        
        # ========== Process each subject ==========
        for i, sub_file in enumerate(subs):
            start_time = time.time()
            print(f"subject {i} : {sub_file}")
            
            # Load timeseries
            timeseries = ft_read_cifti_mod(sub_file)
            
            print(f"Loaded timeseries shape: {timeseries['data'].shape}")
            
            # Check orientation - CIFTI files can be (timepoints, grayordinates) or (grayordinates, timepoints)
            # We need (grayordinates, timepoints) format
            if timeseries['data'].shape[0] < timeseries['data'].shape[1]:
                # Data is (timepoints, grayordinates) - need to transpose
                print(f"Transposing timeseries from {timeseries['data'].shape} to (grayordinates, timepoints)")
                timeseries['data'] = timeseries['data'].T
                # print(f"After transpose: {timeseries['data'].shape}")
            
            # # Check if timeseries and consensus match in size
            # if timeseries['data'].shape[0] != consen_data.shape[0]:
            #     print(f"WARNING: Timeseries grayordinates ({timeseries['data'].shape[0]}) != consensus template ({consen_data.shape[0]})")
            #     if surface_only == 1 and timeseries['data'].shape[0] > consen_data.shape[0]:
            #         print(f"Trimming timeseries to surface only (first {consen_data.shape[0]} grayordinates)")
            #         timeseries['data'] = timeseries['data'][:consen_data.shape[0], :]
                
            
            # ========== Process motion mask ==========
            if power_motion == 1:
                # Load power motion data
                print(f"Loading motion data from power method {B[i]}")
                motion_mat = loadmat(B[i])
                motion_data = motion_mat['motion_data']
                print(f'motion_data shape: {motion_data.shape}')
                orig_motion_filename = os.path.splitext(os.path.basename(B[i]))[0]
                
                # Extract FD thresholds
                allFD = np.zeros(motion_data.shape[1])
                for motion_colum in range(motion_data.shape[1]):
                    fd_val = motion_data[0, motion_colum]['FD_threshold']
                    # print(f"FD threshold: {fd_val}")
                    # Handle nested arrays - keep extracting until we get a scalar
                    while isinstance(fd_val, np.ndarray):
                        if fd_val.size == 1:
                            fd_val = fd_val.item()
                        else:
                            fd_val = fd_val.flatten()[0]
                    allFD[motion_colum] = float(fd_val)
                
                FDidx = np.where(np.round(allFD, 3) == np.round(FD_threshold, 3))[0]
                if len(FDidx) == 0:
                    print(f"Warning: FD threshold {FD_threshold} not found in motion data")
                    print(f"Available FD thresholds: {allFD}")
                    exit(1)
                FDidx = FDidx[0]
                
                # Get frame removal vector
                if combined_outliermask_provided == 1:
                    FDvec_raw = motion_data[0, FDidx]['combined_removal']
                    # print(f"FDvec_raw: {FDvec_raw}")
                else:
                    FDvec_raw = motion_data[0, FDidx]['frame_removal']
                
                # Ensure FDvec is a proper 1D numeric array
                if isinstance(FDvec_raw, np.ndarray):
                    FDvec = np.asarray(FDvec_raw).flatten()
                    # print(f"FDvec after flattening: {FDvec}")
                else:
                    FDvec = np.array(FDvec_raw, dtype=float).flatten()

                
                # Convert any nested arrays to scalars
                if FDvec.dtype == object:
                    if FDvec.size == 1:
                        FDvec = FDvec.item()
                        # print(f"FDvec after cleaning: {FDvec}")
                    else:
                        FDvec = np.concatenate([np.asarray(v).reshape(-1) for v in FDvec])
                        # print(f"FDvec after cleaning: {FDvec}")
                                
                # Ensure FDvec is float to avoid uint8 overflow issues
                FDvec = FDvec.astype(float)
                
                # Convert "1=remove" to "1=keep" (flip 0s and 1s)
                FDvec = 1 - FDvec
                
                # Ensure values are strictly 0 or 1
                FDvec = np.clip(FDvec, 0, 1)
                
                print(f"FDvec after converting: min={FDvec.min()}, max={FDvec.max()}, unique values={np.unique(FDvec)}")
                # print(f"FDvec after converting: {FDvec}")   

                allmasks_before_outliers_removed_FD02.append(FDvec.copy())
                
                # ========== Outlier detection ==========
                if combined_outliermask_provided == 0:
                    if remove_outliers == 1:
                        print("Removal outliers not specified. It will be performed by default.")
                        stdev_temp_filename = f"{file_root_no_ext}_temp.txt"
                        FDvec = CensorBOLDoutliers(wb_command, subs, i, stdev_temp_filename, FDvec)
                    else:
                        print("Motion censoring performed on FD alone. Frames with outliers in BOLD std dev not removed.")
                    
                    good_frames_idx = np.where(FDvec == 1)[0]
                    good_minutes = (len(good_frames_idx) * TR) / 60
                    
                    if good_minutes < 0.5:
                        print(f"Subject {i} has less than 30 seconds of good data")
                        continue
                    elif minutes_to_use > good_minutes:
                        # Use all available frames
                        output_file = os.path.join(project_dir, 
                            f"{orig_motion_filename}_{FD_threshold}_cifti_censor_FD_vector_All_Good_Frames.txt")
                        np.savetxt(output_file, FDvec, fmt='%1.0f')
                    else:
                        # Randomly sample frames to match minutes_to_use
                        good_frames_needed = round(minutes_to_use * 60 / TR)
                        rand_good_frames = sorted(random.sample(list(range(len(good_frames_idx))), good_frames_needed))
                        FDvec_cut = np.zeros(len(FDvec))
                        ones_idx = good_frames_idx[rand_good_frames]
                        FDvec_cut[ones_idx] = 1
                        
                        output_file = os.path.join(project_dir, 
                            f"{orig_motion_filename}_{FD_threshold}_cifti_censor_FD_vector_{minutes_to_use}_minutes_of_data_at_{FD_threshold}_threshold.txt")
                        np.savetxt(output_file, FDvec_cut, fmt='%1.0f')
                        FDvec = FDvec_cut
                else:
                    print("Frames with outliers in BOLD std dev not removed. Maybe you have already performed outlier detection.")
                
                print(f"FDvec: {FDvec.shape}")
                tmask = FDvec.astype(bool)
                tmask = np.asarray(tmask).reshape(-1)
                print(f"tmask: {tmask.shape}")
            else:
                # Load motion mask from file
                if isinstance(B[i], (np.ndarray, list)):
                    tmask = np.array(B[i]).astype(bool)
                    FDvec = tmask.astype(float)
                else:
                    tmask = np.loadtxt(B[i]).astype(bool)
                    FDvec = tmask.astype(float)
            
            allmasks_outliers_removed_FD02.append(FDvec)
            
            # ========== Censor timeseries ==========
            timeseries_data = timeseries['data']
            
            # Debug: print shapes
            print(f"Timeseries shape: {timeseries_data.shape}")
            print(f"tmask shape: {tmask.shape}")
            print(f"tmask length: {len(tmask)}")
            print(f"Number of True values in tmask: {np.sum(tmask)}")
            
            ## Need to ensure if this check is actually working
            # Ensure tmask length matches number of timepoints
            if len(tmask) != timeseries_data.shape[1]:
                print(f"WARNING: tmask length ({len(tmask)}) does not match number of timepoints ({timeseries_data.shape[1]})")
                if len(tmask) < timeseries_data.shape[1]:
                    print(f"ERROR: tmask is too short. Skipping subject {i}.")
                    subs.pop(i)
                    B.pop(i)
                    print(f"  After popping: len(subs) = {len(subs)}, len(B) = {len(B)}")
                    continue
                elif len(tmask) > timeseries_data.shape[1]:
                    print(f"Truncating tmask to match timepoints")
                    tmask = tmask[:timeseries_data.shape[1]]
            
            # Apply censoring
            timeseries_data = timeseries_data[:, tmask]
            print(f"Censored timeseries shape: {timeseries_data.shape}")
            
            # ========== Compute correlations for each network ==========
            corrs = np.zeros((timeseries['data'].shape[0], len(network_names)))
            
            for j, net_name in enumerate(network_names):
                if j in [3, 5, 16]:  # Skip empty network slots
                    continue
                
                # Get network mask (MATLAB is 1-indexed)
                inds = consen_data == (j + 1)
                net_idx = np.where(inds)[0]
                
                # Debug
                print(f"Processing network {j} ({net_name})")
                print(f"  Network mask shape: {net_idx.shape}, sum: {np.sum(inds)}")
                print(f"  Timeseries data shape: {timeseries_data.shape}")
                
                # Calculate network average - inds indexes grayordinates (rows), not timepoints
                # timeseries_data is (grayordinates, timepoints)
                # inds is (grayordinates,)
                subNetAvg = np.nanmean(timeseries_data[net_idx, :], axis=0)
                print(f"  Network average shape: {subNetAvg.shape}")
                
                # Compute correlations
                if surface_only == 1:
                    voxel_range = range(ncortgrey)
                else:
                    voxel_range = range(len(timeseries_data))
                
                for voxel in voxel_range:
                    goodvox = ~np.isnan(timeseries_data[voxel, :])
                    if np.sum(goodvox) > 1:
                        corrs[voxel, j] = paircorr_mod(
                            subNetAvg[goodvox], 
                            timeseries_data[voxel, goodvox]
                        )
            
            seedmapstimeseries.append(corrs)
            elapsed = time.time() - start_time
            print(f"Elapsed time: {elapsed:.2f} seconds")
        
        # ========== Save seedmaps ==========
        print(f"\nSaving seedmaps for {len(seedmapstimeseries)} subjects...")
        # Convert to cell array format for MATLAB compatibility
        seedmaps_cell = np.empty((len(seedmapstimeseries), 1), dtype=object)
        for i, sm in enumerate(seedmapstimeseries):
            seedmaps_cell[i, 0] = sm
        
        masks_cell = np.empty((len(allmasks_outliers_removed_FD02), 1), dtype=object)
        for i, mask in enumerate(allmasks_outliers_removed_FD02):
            masks_cell[i, 0] = mask
        
        masks_before_cell = np.empty((len(allmasks_before_outliers_removed_FD02), 1), dtype=object)
        for i, mask in enumerate(allmasks_before_outliers_removed_FD02):
            masks_before_cell[i, 0] = mask
        
        savemat(
            seedmaps_file,
            {
                'seedmapstimeseries': seedmaps_cell,
                'allmasks_outliers_removed_FD02': masks_cell,
                'allmasks_before_outliers_removed_FD02': masks_before_cell
            }
        )
    
    # ========== Check for NaNs ==========
    print("Checking for nans in seeds...")
    print(f"Number of subjects to check: {len(subs)}")
    print(f"Number of seedmap subjects: {len(seedmapstimeseries)}")
    badsubidx = []
    cleansubs = subs.copy()
    
    for i in range(len(subs)):
        has_nan = False
        for j, net_name in enumerate(network_names):
            if j in [3, 5, 16]:  # Skip empty slots
                continue
            
            if np.sum(np.isnan(seedmapstimeseries[i][:, j])) > 0:
                print(f"This subject {i} has Nans")
                print(f"File with nans is: {subs[i]}")
                print(f"Network with nans is: {net_name}")
                print(f"Number of greyordinates with nans = {np.sum(np.isnan(seedmapstimeseries[i][:, j]))}")
                badsubidx.append(i)
                has_nan = True
                break
    
    if len(badsubidx) == 0:
        print("Congratulations, your data has no nans.")
    else:
        badsubidx = sorted(list(set(badsubidx)))
        cleansubs = [s for idx, s in enumerate(subs) if idx not in badsubidx]
        seedmapstimeseries = [s for idx, s in enumerate(seedmapstimeseries) if idx not in badsubidx]
        print(f"{len(badsubidx)} subjects were removed from average for 'NaNs' in greyordinates.")
    
    # ========== Average across subjects ==========
    avgSeedmaps = [None] * len(network_names)
    seed_matrix = np.zeros((seedmapstimeseries[0].shape[0], len(network_names)))
    
    for j, net_name in enumerate(network_names):
        print(f"Processing network {j} ({net_name})")
        if j in [3, 5, 16]:  # Skip empty slots
            print(f"Skipping empty network {j}")
            continue
        
        # Fisher transform, average, and inverse transform
        grpNetAve = np.zeros(seedmapstimeseries[0].shape[0])
        for i in range(len(cleansubs)):
            grpNetAve += np.arctanh(seedmapstimeseries[i][:, j])
        
        grpNetAve /= len(cleansubs)
        avgSeedmap = np.tanh(grpNetAve)
        
        # ========== Z-score within regions ==========
        if Zscore_regions == 1:
            print("Converting to Zscores")
            try:
                avgSeedmap[0:29696] = stats.zscore(avgSeedmap[0:29696])
                avgSeedmap[29697:59412] = stats.zscore(avgSeedmap[29697:59412])
                if surface_only == 0:
                    avgSeedmap[59413:91282] = stats.zscore(avgSeedmap[59413:91282])
            except Exception as e:
                print(f"Error z-scoring network {j}: {e}")
                return
        
        avgSeedmaps[j] = avgSeedmap
        seed_matrix[:, j] = avgSeedmap
        
        # ========== Save individual network map ==========
        temp_file_path = cleansubs[i]
        temp_img_data = ciftiopen(temp_file_path, wb_command)
        print(f"Loaded temp_img: {temp_img_data}")
        # Use 'img' key to get the full Cifti2Image object
        if 'img' in temp_img_data:
            temp_cifti = temp_img_data['img']
        else:
            temp_cifti = temp_img_data['cdata']

        print(f"Original temp_img shape: {temp_cifti.shape}")
        
        if surface_only == 1:
            template_path = settings['path'][10]  # Adjust index for surface-only template
            temp_img_data = ciftiopen(template_path, wb_command)
            if 'img' in temp_img_data:
                temp_cifti = temp_img_data['img']
            else:
                temp_cifti = temp_img_data['cdata']
            temp_cifti = np.zeros((temp_cifti.shape[0], 1))
            print(f"Template shape: {temp_cifti.shape}")
        
        # Get original data to understand shape
        original_data = temp_cifti
        print(f"Original data shape: {original_data.shape}")

        # Create new data with correct shape (1, 91282)
        new_data = avgSeedmap.reshape(1, -1)
        print(f"New data shape: {new_data.shape}")

        # Get the spatial axis (brain models axis) from template
        brain_models_axis = temp_cifti.header.get_axis(1)

        # Create new scalar axis for single map
        scalar_axis = nib.cifti2.ScalarAxis([net_name])

        # Create new header with correct dimensions
        new_header = nib.cifti2.Cifti2Header.from_axes((scalar_axis, brain_models_axis))

        # Create new CIFTI image with your data
        new_cifti = nib.Cifti2Image(new_data, new_header)
        print(f"New CIFTI shape: {new_cifti.shape}")
        print(f"New CIFTI : {new_cifti}")
         
        # Save CIFTI file
        if Zscore_regions == 1:
            output_path = os.path.join(project_dir, 
                f'seedmaps_{file_root_no_ext}_{net_name}_networkZscored.dtseries.nii')
        else:
            output_path = os.path.join(project_dir, 
                f'seedmaps_{file_root_no_ext}_{net_name}_network.dtseries.nii')
        
        ciftisave(new_cifti, output_path)
    
    # ========== Save all maps ==========
    num_cleaned_subs = len(cleansubs)
    print(f'cleaned subs: {type(cleansubs)}')

    # Convert lists to MATLAB cell array format (1 x N)
    cleansubs_cell = np.empty((1, len(cleansubs)), dtype=object)
    cleansubs_cell[0, :] = cleansubs

    B_cell = np.empty((1, len(B)), dtype=object)
    B_cell[0, :] = B

    subs_cell = np.empty((1, len(subs)), dtype=object)
    subs_cell[0, :] = subs
    
    # Convert indices to double (float64) for MATLAB compatibility
    bad_subs_idx_double = np.array(bad_subs_idx, dtype=np.float64) if 'bad_subs_idx' in locals() else np.array([], dtype=np.float64)
    good_subs_idx_double = np.array(good_subs_idx, dtype=np.float64) if 'good_subs_idx' in locals() else np.array([], dtype=np.float64)


    save_params = {
        'B': B_cell,  # Convert to cell array
        'seed_matrix': seed_matrix,
        'subs': subs_cell,  # Convert to cell array
        'Zscore_regions': Zscore_regions,
        'power_motion': power_motion,
        'remove_outliers': remove_outliers,
        'surface_only': surface_only,
        'use_only_subjects_that_pass_motion_criteria': use_only_subjects_that_pass_motion_criteria,
        'combined_outliermask_provided': combined_outliermask_provided,
        'include_scan_net': include_scan_net,
        'bad_subs_idx': bad_subs_idx_double if 'bad_subs_idx' in locals() else [],
        'good_subs_idx': good_subs_idx_double if 'good_subs_idx' in locals() else [],
        'cleansubs': cleansubs_cell,  # Convert to cell array
        'FD_threshold': FD_threshold,
        'minutes_to_use': minutes_to_use
    }

    
    if Zscore_regions == 1:
        output_file = os.path.join(project_dir, 
            f'seedmaps_{file_root_no_ext}_n_{num_cleaned_subs}_all_networksZscored.mat')
    else:
        output_file = os.path.join(project_dir, 
            f'seedmaps_{file_root_no_ext}_n_{num_cleaned_subs}_all_networks.mat')
    
    savemat(output_file, save_params)
    
    print("Done making network templates based on the subjects you provided.")


# ============================================================================
# Helper Functions
# ============================================================================

def settings_comparematrices():
    """
    Load settings with paths and configuration for CIFTI template processing.
    
    Returns:
    --------
    dict : Dictionary with paths and configuration
    """
    settings = {
        'path': [
            '/projects/standard/faird/shared/code/external/utilities/gifti-1.6',  # 0 (MATLAB: path{1})
            '/projects/standard/faird/shared/code/internal/utilities/Matlab_CIFTI',  # 1 (MATLAB: path{2})
            '/projects/standard/faird/shared/code/external/utilities/Matlab_effect_size_toolbox/',  # 2 (MATLAB: path{3})
            '/projects/standard/faird/shared/code/internal/utilities/hcp_comm_det_damien/Gordan_subcortical_template_ptseries.conc_AVG.pconn.nii',  # 3 (MATLAB: path{4})
            '/projects/standard/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/Networks_template_cleaned.pscalar.nii',  # 4 (MATLAB: path{5})
            '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/Networks_template_cleaned.dscalar.nii',  # 5 (MATLAB: path{6})
            '/projects/standard/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/120_LR_minsize400_recolored_manualconsensus4.dconn.nii',  # 6 (MATLAB: path{7})
            '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dscalar.nii',  # 7 (MATLAB: path{8})
            '/projects/standard/faird/shared/code/internal/utilities/hcp_comm_det_damien/Merged_HCP_best80_dtseries.conc_AVG.dconn.nii',  # 8 (MATLAB: path{9})
            '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates.dtseries.nii',  # 9 (MATLAB: path{10})
            '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dtseries.nii',  # 10 (MATLAB: path{11})
            '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/91282_Greyordinates_surf_only.dscalar.nii',  # 11 (MATLAB: path{12})
            '/projects/standard/faird/shared/code/internal/utilities/figure_maker/MSC01_template_quad_scaled_v3_legend_fixed_MSI.scene',  # 12 (MATLAB: path{13})
            '/projects/standard/faird/shared/code/internal/utilities/figure_maker/MSC01_template_scene_subcort_scalar_MSI.scene',  # 13 (MATLAB: path{14})
            '/projects/standard/faird/shared/code/internal/utilities/figure_maker/quick_network_pic.sh',  # 14 (MATLAB: path{15})
            '/projects/standard/faird/shared/code/internal/utilities/community_detection/fair/supporting_files/EUGEODistancematrix_XYZ_255interhem_unit8.mat',  # 15 (MATLAB: path{16})
        ],
        'path_wb_c': '/projects/standard/feczk001/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command',
        'path_template_nets': '/projects/standard/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/support_files/seedmaps_ADHD_smoothed_dtseries315_all_networks_Zscored.mat'
    }
    
    return settings


def ft_read_cifti_mod(filepath):
    """
    Read CIFTI file using nibabel.
    
    Parameters:
    -----------
    filepath : str
        Path to CIFTI file
        
    Returns:
    --------
    dict : Dictionary with 'data' and 'header'
    """
    img = nib.load(filepath)
    return {
        'data': img.get_fdata().squeeze(),
        'header': img.header,
        'img': img
    }


def check_twins_motion(motion_file, TR, project_dir, output_name, 
                       FD_column, FD_conc_file, get_mean_FD, 
                       use_outlier_mask, split_motion, splits):
    """
    Check motion data and return minutes of good data.
    
    Parameters:
    -----------
    motion_file : str or list
        Path to motion file (.mat or .txt)
    TR : float
        Repetition time in seconds
    project_dir : str
        Output directory
    output_name : str
        Output filename
    FD_column : int
        Which FD column to use (e.g., 21 for FD=0.2)
    FD_conc_file : str or None
        Path to FD concatenated file
    get_mean_FD : int
        Whether to get mean FD (0 or 1)
    use_outlier_mask : int
        Use outlier detection mask if available (0 or 1)
    split_motion : int
        Whether to split motion (0 or 1)
    splits : int
        Number of splits
        
    Returns:
    --------
    tuple : (all_subjects_minutes, all_mean_FD)
    """
    
    # Handle single file or list
    if isinstance(motion_file, str):
        if motion_file.endswith('.conc'):
            with open(motion_file, 'r') as f:
                all_motion_conc = [line.strip() for line in f.readlines() if line.strip()]
        else:
            all_motion_conc = [motion_file]
    elif isinstance(motion_file, list):
        all_motion_conc = motion_file
    else:
        all_motion_conc = [motion_file]
    
    # Check file extension
    _, ext = os.path.splitext(all_motion_conc[0])
    
    # Initialize arrays
    good_frames = np.zeros(len(all_motion_conc))
    all_mean_FD = np.zeros(len(all_motion_conc))
    all_subjects_minutes = np.zeros(len(all_motion_conc))
    all_possible_minutes = np.zeros(len(all_motion_conc))
    all_poss_frames = np.zeros(len(all_motion_conc))
    
    # Process each subject
    for i in range(len(all_motion_conc)):
        # print(f"Loading subject {i + 1}")
        
        if ext == '.txt':
            # Load text file with 0s and 1s
            FD_vector_1isbad = np.loadtxt(all_motion_conc[i])
            good_frames[i] = np.sum(FD_vector_1isbad)
            all_mean_FD[i] = 0
            
            # Calculate minutes
            possible_frames = len(FD_vector_1isbad)
            good_frames_in_minutes = good_frames[i] * TR / 60
            possible_minutes = possible_frames * TR / 60
            
        else:
            # Load .mat file
            try:
                motion_mat = loadmat(all_motion_conc[i])
                motion_data = motion_mat['motion_data']
                
                # Get frame removal vector
                if use_outlier_mask == 1:
                    try:
                        # Try to use combined_removal mask
                        # print(f"Subject {i} has outlier detection mask. Using combined_removal.")
                        # print(f"  Subject {i}: motion_data shape={motion_data.shape}")
                        FD02_vector = motion_data[0, FD_column - 1]['combined_removal'].flatten()[0]
                        # print(f"  Subject {i}: combined_removal length={FD02_vector}")
                    except (KeyError, IndexError):
                        print(f"Subject {i} does not have outlier detection mask. Using frame_removal instead.")
                        FD02_vector = motion_data[0, FD_column - 1]['frame_removal'].flatten()
                else:
                    FD02_vector = motion_data[0, FD_column - 1]['frame_removal'].flatten()
                
                # Handle nested arrays
                if FD02_vector.dtype == object:
                    FD02_clean = np.zeros(FD02_vector.shape[0])
                    for idx, val in enumerate(FD02_vector):
                        if isinstance(val, np.ndarray):
                            FD02_clean[idx] = val.item() if val.size == 1 else val.flatten()[0]
                        else:
                            FD02_clean[idx] = float(val)
                    FD02_vector = FD02_clean
                
                # MATLAB: frame_removal has 1=remove, 0=keep
                # We need to count frames to KEEP (good frames)
                # So good_frames = sum(abs(FD02_vector - 1))
                # This converts 1 -> 0, and 0 -> 1, then sums the 1s (frames to keep)
                # Ensure FDvec is float to avoid uint8 overflow issues
                FD02_vector = FD02_vector.astype(float)
                
                # Convert "1=remove" to "1=keep" (flip 0s and 1s)
                FD02_vector = 1 - FD02_vector
                
                # Ensure values are strictly 0 or 1
                FD02_vector = np.clip(FD02_vector, 0, 1)
                good_frames[i] = np.sum(FD02_vector)
                # print(f"  Subject {i}: FD vector={FD02_vector}")
                
                print(f"  Subject {i+1}: FD vector length={len(FD02_vector)}, good frames={good_frames[i]}")
                
                # Get mean FD if requested
                if get_mean_FD == 1:
                    try:
                        all_mean_FD[i] = float(motion_data[0, FD_column - 1]['remaining_frame_mean_FD'])
                    except (KeyError, IndexError):
                        all_mean_FD[i] = 0
                
                # Get TR from motion data or use provided TR
                try:
                    epi_TR = float(motion_data[0, FD_column - 1]['epi_TR'])
                except (KeyError, IndexError, ValueError):
                    epi_TR = TR
                
                # Calculate minutes
                possible_frames = len(FD02_vector)
                good_frames_in_minutes = good_frames[i] * epi_TR / 60
                possible_minutes = possible_frames * epi_TR / 60
                
            except Exception as e:
                print(f"Error loading motion file {i}: {e}")
                good_frames[i] = 0
                all_mean_FD[i] = 0
                good_frames_in_minutes = 0
                possible_minutes = 0
                possible_frames = 0
        
        # Store results
        all_subjects_minutes[i] = good_frames_in_minutes
        all_possible_minutes[i] = possible_minutes
        all_poss_frames[i] = possible_frames
        
        print(f"Subject {i+1}: {good_frames_in_minutes:.2f} minutes available out of {possible_minutes:.2f} possible minutes")
    
    # Save summary
    summary_data = {
        'good_frames': good_frames,
        'all_subjects_minutes': all_subjects_minutes,
        'all_mean_FD': all_mean_FD,
        'all_poss_frames': all_poss_frames,
        'all_possible_minutes': all_possible_minutes,
        'all_motion_conc_file': motion_file,
        'all_motion_conc': all_motion_conc
    }
    
    summary_file = os.path.join(project_dir, f'{output_name}_summary.mat')
    savemat(summary_file, summary_data)
    
    print('Done checking motion')
    
    return all_subjects_minutes, all_mean_FD


def CensorBOLDoutliers(wb_command, subs, idx, stdev_temp_filename, FDvec):
    """
    Detect and censor BOLD outliers using standard deviation method.
    
    This is a placeholder - implement based on your outlier detection logic.
    
    Parameters:
    -----------
    wb_command : str
        Path to wb_command
    subs : list
        List of subject files
    idx : int
        Current subject index
    stdev_temp_filename : str
        Temporary filename for output
    FDvec : ndarray
        Frame removal vector
        
    Returns:
    --------
    ndarray : Updated frame removal vector
    """
    # Placeholder implementation
    # TODO: Implement actual outlier detection
    print(f"Performing outlier detection for subject {idx}")
    return FDvec


def paircorr_mod(x, y):
    """
    Calculate Pearson correlation between two vectors.
    
    Parameters:
    -----------
    x, y : ndarray
        Input vectors
        
    Returns:
    --------
    float : Correlation coefficient
    """
    if len(x) < 2 or len(y) < 2:
        return np.nan
    
    # Remove any remaining NaNs
    valid_idx = ~(np.isnan(x) | np.isnan(y))
    if np.sum(valid_idx) < 2:
        return np.nan
    
    return np.corrcoef(x[valid_idx], y[valid_idx])[0, 1]


def ciftiopen(filepath, wb_command):
    """
    Open CIFTI file - mimics MATLAB's ciftiopen.
    
    Parameters:
    -----------
    filepath : str
        Path to CIFTI file
    wb_command : str
        Path to wb_command (not used in Python version)
        
    Returns:
    --------
    dict : Dictionary with CIFTI data and metadata
    """
    img = nib.load(filepath)
    return {
        'cdata': img.get_fdata(),
        'header': img.header,
        'img': img
    }


def ciftisave(cifti_struct, filepath):
    """
    Save CIFTI file - mimics MATLAB's ciftisave.
    
    Parameters:
    -----------
    cifti_struct : dict or Cifti2Image
        Dictionary with CIFTI data (with 'cdata', 'header', 'img' keys) or Cifti2Image object
    filepath : str
        Output filepath
    wb_command : str, optional
        Path to wb_command (not used in Python version, kept for compatibility)
    """
    # Handle different input types
    if isinstance(cifti_struct, nib.Cifti2Image):
        # Direct Cifti2Image object
        new_img = cifti_struct
    elif isinstance(cifti_struct, dict):
        if 'img' in cifti_struct and isinstance(cifti_struct['img'], nib.Cifti2Image):
            # Use existing Cifti2Image from dict
            new_img = cifti_struct['img']
        elif 'cdata' in cifti_struct:
            # Need to reconstruct - check if we have a template image
            if 'img' in cifti_struct:
                # Reconstruct with updated data
                img = cifti_struct['img']
                new_img = nib.Cifti2Image(cifti_struct['cdata'], header=img.header)
            elif 'header' in cifti_struct:
                # Use provided header
                new_img = nib.Cifti2Image(cifti_struct['cdata'], header=cifti_struct['header'])
            else:
                raise ValueError("cifti_struct must contain 'img' or 'header' key")
        else:
            raise ValueError("cifti_struct dict must contain 'cdata' or 'img' key")
    else:
        raise TypeError(f"cifti_struct must be dict or Cifti2Image, got {type(cifti_struct)}")
    
    # Save the image
    nib.save(new_img, filepath)
    print(f"Saved CIFTI file to: {filepath}")


# ============================================================================
# Main execution
# ============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Create CIFTI templates with motion censoring')
    parser.add_argument('--timeseries', required=True, help='Path to timeseries conc file')
    parser.add_argument('--TR', type=float, required=True, help='Repetition time in seconds')
    parser.add_argument('--motion', required=True, help='Path to motion conc file')
    parser.add_argument('--project_dir', required=True, help='Output directory')
    parser.add_argument('--zscore_regions', type=int, default=0, help='Z-score within regions (0 or 1)')
    parser.add_argument('--power_motion', type=int, default=1, help='Use Power motion method (0 or 1)')
    parser.add_argument('--remove_outliers', type=int, default=1, help='Remove outlier frames (0 or 1)')
    parser.add_argument('--surface_only', type=int, default=0, help='Surface only (0 or 1)')
    parser.add_argument('--use_motion_criteria', type=int, default=1, help='Filter by motion (0 or 1)')
    parser.add_argument('--combined_mask', type=int, default=0, help='Combined mask provided (0 or 1)')
    parser.add_argument('--include_scan', type=int, default=0, help='Include SCAN network (0 or 1)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.project_dir, exist_ok=True)
    
    # Run the main function
    makeCiftiTemplates_RH(
        dt_or_ptseries_conc_file=args.timeseries,
        TR=args.TR,
        all_motion_conc_file=args.motion,
        project_dir=args.project_dir,
        Zscore_regions=args.zscore_regions,
        power_motion=args.power_motion,
        remove_outliers=args.remove_outliers,
        surface_only=args.surface_only,
        use_only_subjects_that_pass_motion_criteria=args.use_motion_criteria,
        combined_outliermask_provided=args.combined_mask,
        include_scan_net=args.include_scan
    )
