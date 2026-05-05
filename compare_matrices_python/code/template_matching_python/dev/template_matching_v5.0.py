import os
import subprocess
import sys
import nibabel as nb
from nibabel import cifti2
import scipy.io as sio       # scipy is for accessing mat files
import numpy as np
import time
import pandas as pd

# TO DO:
# i. warning if input ciftis and templates are not transformed correctly (both pearson, fisher, or z-scored)

# update log
# v3.0 
# i.   updated eta_to_vox_modified = eta_to_template_vox.copy() to prevent original from overwrite
# ii.  rather than returning NaN if cmap == 0 or tmap == 0 now only returns NaN if SStot = 0
# iii. tried new method of writing dscalar file
# v4.0 
# major updates:
# i.  returned to returning NaN if cmap == 0 or tmap == 0 as this matches matlab code functionality
#     confirmed matching functionality by modifying matlab code to save both eta matrices
# ii. debugged method of writing new cifti file to properly save out new network assignments
# minor updates:
# i.  transforming network_names with np.array(network_names,dtype=object) 
#     now saving network assignments similarly in the .mat file
# v5.0.
# major updates:
# i.  saving out both dscalars (with and without network refinement)
# ii. setting up to run with an sbatch command
# minor updates:
# i.   cleaned up and deleted old code
# ii.  user warning if dconn and template have different number of greyordinates
# iii. updated the name of refine_scan_network step to refine_template_matching
#      changed the name of refinement steps to reflect that refinement happens for multiple networks (not just scan)
#      includes: SMd, SMl, SCAN, DMN, PMN, PON, SAL, CO
# iv.  refine_template_matching is now a user-specified shell script option

# get the current time
totaltimer = time.time()

# get arguements from command line
# input dconn
dconn_infile = [sys.argv[1]][0]
# input template and refinement options
template_infile = [sys.argv[2]][0]
template_minthreshold = [sys.argv[3]][0]
refine_template_matching = [sys.argv[4]][0]
refinement_minthreshold = [sys.argv[5]][0]
# output matlab file
matfile_out = [sys.argv[6]][0]
# output dscalar files
dscalar_template_infile = [sys.argv[7]][0]
dscalar_outfile = [sys.argv[8]][0]
dscalar_outfile_refined = [sys.argv[9]][0]

# user specified options
TEMPLATEMINIMUM = float(template_minthreshold)
REFINEMENTTEMPLATEMINIMUM = float(refinement_minthreshold)
network_names = ['DMN','Vis','FP','','DAN','','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON','','SCAN']

# can hardcode paths to necessary directories, files and code
# this makes it easier to debug the code and run on command line:

# wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/'
# dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn'
# cifti_output_folder='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching'

# input files
# dconn_infile = '/home/btervocl/shared/projects/MINT/kgodfrey/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH.dconn.nii'
# template_infile='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/support_files/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat'
# dscalar_template_infile='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/support_files/sub-10227_ses-combined_task-cross_grims_recoloredCifti2.dscalar.nii'

# output files
# matfile_out='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/v6_interpolated/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.mat'
# dscalar_outfile='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/v6_interpolated/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.dscalar.nii'
# dscalar_outfile_refined='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/v6_interpolated/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonScanTM.dscalar.nii'

# options
# TEMPLATEMINIMUM = float(1.00)
# REFINEMENTTEMPLATEMINIMUM = float(3.00)
# network_names = ['DMN','Vis','FP','','DAN','','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON','','SCAN']

file_checker=True
if file_checker:
    # check for input dconn
    if os.path.isfile(dconn_infile) == False:
        print('your dconn path doesnt exist')
        print('\tinput dconn file path provided: ' + dconn_infile)
    # check for input template
    if os.path.isfile(template_infile) == False:
        print('your template path doesnt exist')
        print('\tinput template file path provided: ' + template_infile)
    # check for dscalar template
    if os.path.isfile(template_infile) == False:
        print('your dscalar template path doesnt exist')
        print('\tinput dscalar template file path provided: ' + dscalar_template_infile)
        print('\tscript cant produce output .dscalar.nii')

load_data=True
if load_data:
    print('loading data')
    print('\t')
    steptimer = time.time()
    
    print('\tloading the correlation matrix')
    subject_cii = nb.load(dconn_infile)
    corr_mat_full=subject_cii.get_fdata()  # get dconn out of input cifti
    greyordinates = corr_mat_full.shape[0]
    del subject_cii # save memory
    print('\tnumber of greyordinates: ' + str(greyordinates))
    print('\t ')
    
    print('\tloading the template')
    cifti_template = sio.loadmat(template_infile)
    cifti_template_mat_full = cifti_template.get("seed_matrix")
    cifti_template_mat_full = cifti_template_mat_full.copy()
    cifti_template_mat_full[cifti_template_mat_full <= TEMPLATEMINIMUM] = 'NaN'
    del cifti_template
    print('\tnumber of template greyordinates: ' + str(cifti_template_mat_full.shape[0]))
    print('\ttemplate minimum threshold set at: ' + str(TEMPLATEMINIMUM))
    print('\ttemplate greyordinates below minimum set to NaN')
    print('\t')
    
    if (corr_mat_full.shape[0] != cifti_template_mat_full.shape[0]):
        print('\t WARNING: INPUT DCONN AND INPUT TEMPLATE DIFFER IN NUMBER OF GREYORDINATES')
        print('\t unequal number of greyordinates will cause errors with template matching')
        print('\t check resolution of inputs or inclusion/exclusion of subcortical greyordinates')
        print('\t')
    
    print('\tloading data complete')
    print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
    print('\t ')

run_template_matching=True
if run_template_matching:
    
    print('template matching at threshold: ' + str(TEMPLATEMINIMUM))
    print('\t')
    steptimer = time.time()
    
    # create empty array to store template matched outputs (1 row per greyordinate and 1 column per network)
    # eta_to_template_vox = calculate association to template using RH sum of squares
    # r_to_template_vox = calculate association to template with Pearson correlation
    eta_to_template_vox = np.zeros((int(len(corr_mat_full)),int(len(network_names))))
    r_to_template_vox = np.zeros((int(len(corr_mat_full)),int(len(network_names))))
    
    print('\tcalculating similarity of greyordinates to template')
    print('\t')
    
    # i will loop through each vertex from the correlation matrix
    for i in range(greyordinates):
        
        if i % 5000 == 0: # print a progress update every 5000 voxels
            print('\tcalculating greyordinate: ' + str(i))
            print('\ttemplate matching has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
        for j in range(len(network_names)): # j will loop through each of the candidate networks to match to
            if j in (3,5,16):
                pass
            else:
                # identify greyordinates which are not 'NaN' in cmap and tmap
                # 'NaN' greyordinates specified earlier w/ user-specified tmap threshold
                goodvox = np.where(~np.isnan(corr_mat_full[i,:]) & ~np.isnan(cifti_template_mat_full[:,j]))
                # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                cmap = corr_mat_full[i,goodvox]
                # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                tmap = cifti_template_mat_full[goodvox,j]
                
                # calculation of eta values using RH within and between sum of squares ##########################################
                # calculate the grand mean across all elements of tmap and cmap
                Mgrand  = ((np.sum(tmap)/(tmap.shape[1])) + (np.sum(cmap)/(cmap.shape[1]))) / 2
                # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                Mwithin = np.mean([tmap, cmap], axis=0)
                # get a sum of squares within
                SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                # get a sum of squares total
                SStot    = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                # eta is calculated as 1 - this proportional variance explained
                eta = (1 - (SSwithin/SStot))
                eta_to_template_vox[i,j] = (1 - (SSwithin/SStot))
                
                # calculation of r values using correlation ######################################################################
                r_to_template_vox[i,j] = np.corrcoef(cmap, tmap)[0,1]
                # delete outputs before the next loop
                del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox
    
    # winner-take-all assignment of greyordinates to networks
    print('\t')
    print('\twinner-take-all assignment of greyordinates to networks')
    print('\t')
    # get new network assignments by using idxmax to return which column has the highest eta
    new_subject_labels = np.nanargmax(eta_to_template_vox,axis=1)
    # get new network assignments by using idxmax to return which column has the highest r 
    new_subject_labels_r = np.nanargmax(r_to_template_vox,axis=1)
    # wrap up this step
    print('\ttemplate matching complete')
    print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
    print('\t')

# save template matching results from this first step
# hardcoded, but might be useful for debugging
save_first_step_mat=False
if save_first_step_mat:
    mat_dict={"network_names": network_names
            ,"eta_to_template_vox": eta_to_template_vox
            ,"new_subject_labels":new_subject_labels}
    matfile_out='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matchingv2.0_Z/sub-10227_ses-combined_task-cross_noFisher_zscoreRH2.mat'    
    sio.savemat(matfile_out,mat_dict,appendmat=False)

# get eta_to_template_vox from a previously created .mat file 
# hardcoded, but might be useful for debugging or to load priors
load_old_mat = False
if load_old_mat:
    print('loading old mat file')
    mat_filepath='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matchingv2.0_Z/sub-10227_ses-combined_task-cross_TM.mat'
    mat_data = sio.loadmat(mat_filepath)
    eta_to_template_vox = np.array(mat_data.get('eta_to_template_vox'))
    eta_to_template_vox_modified = np.array(mat_data.get('eta_to_template_vox_modified'))
    new_subject_labels = np.array(mat_data.get('new_subject_labels'))
    new_subject_labels_refined = np.array(mat_data.get('new_subject_labels_refined'))

# optional refinement of template matching with new threshold
# currently specified at command line
if refine_template_matching:
    print('WARNING: refinement not currently recommended')
    print(' ')
    print('template matching refinement at threshold: ' + str(REFINEMENTTEMPLATEMINIMUM))
    print('')
    refinetimer = time.time()
    
    # copy eta results into a new array
    eta_to_template_vox_modified = eta_to_template_vox.copy()
    r_to_template_vox_modified = r_to_template_vox.copy()
    
    # copy winner-take-all results into new array
    new_subject_labels_refined=new_subject_labels.copy()
    new_subject_labels_r_refined=new_subject_labels_r.copy()
    
    # get the seed matrix from before but apply the new threshold
    cifti_template = sio.loadmat(template_infile)
    cifti_template = cifti_template.get("seed_matrix")
    cifti_template_mat_full_refined = cifti_template.copy()
    #cifti_template_mat_full_refined = cifti_template.get("seed_matrix")
    cifti_template_mat_full_refined[cifti_template_mat_full_refined <= REFINEMENTTEMPLATEMINIMUM] = 'NaN'
    
    # refine greyordinates in SMd(10), SMl(11), or SCAN (18) networks
    refine_motor_grays=True
    if refine_motor_grays:
        
        print('\trefinement of motor grays:')
        steptimer=time.time()
        
        # find greyordinates belonging to SMd(10), SMl(11), or SCAN (18) but index w/ 9,10,17 b/c python is base 0
        motor_grays=np.transpose(np.where((new_subject_labels_refined == 9) | (new_subject_labels_refined == 10) | (new_subject_labels_refined == 17)))
        
        # run template matching looping over all motor_grays and all networks
        for i in range(motor_grays.shape[0]):
            
            if i % 5000 == 0:
                print('\t\tmotor grey refinement has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
                
            for j in range(eta_to_template_vox_modified.shape[1]):
                
                # skip networks that don't exist
                if j in (3,5,16):
                    pass
                
                else:
                    
                    # identify voxels which are not 'NaN' in cmap and tmap
                    motor_gray = motor_grays[i].item()
                    goodvox = np.where(~np.isnan(corr_mat_full[motor_gray,:]) & ~np.isnan(cifti_template_mat_full_refined[:,j]))
                    # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                    cmap = corr_mat_full[motor_gray,(goodvox)]
                    # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                    tmap = cifti_template_mat_full_refined[goodvox,j]
                    
                    # check if cmap or tmap are zero, and if so return 'nan'
                    if (np.sum(cmap) == 0) and (np.sum(tmap) == 0):
                        eta_to_template_vox_modified[motor_gray,j] = float("nan")
                        r_to_template_vox_modified[i,j] = float("nan")
                    
                    # if cmap and tmap are not zero calculate effect size as usual
                    else:
                        
                        # calculation of eta values using RH within and between sum of squares ##########################################
                        # calculate the grand mean across all elements of tmap and cmap
                        Mgrand = (np.mean(tmap) + np.mean(cmap)) / 2
                        # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                        Mwithin = np.mean([tmap, cmap], axis=0)
                        # get a sum of squares within
                        SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                        # get a sum of squares total
                        SStot    = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                        # eta is calculated as 1 - this proportional variance explained
                        eta = (1 - (SSwithin/SStot))
                        eta_to_template_vox_modified[motor_gray,j] = (1 - (SSwithin/SStot))
                        
                        # calculation of r values using correlation ######################################################################
                        r_to_template_vox_modified[motor_gray,j] = np.ma.corrcoef(cmap, tmap)[0,1]
                        # delete outputs before the next loop
                        del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox
        
        print('\t\trefinement of motor grays complete')
        print('\t\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
        print('\t\t')
    
    # refine greyordinates in DMN(1), PMN(15), or PON(16)
    refine_default_grays=False
    if refine_default_grays:
        
        print('\trefinement of default mode grays:')
        steptimer=time.time()
        
        # find the networks that belong to DMN(1), PMN(15), or PON (16), using 0,14,15 b/c base 0
        dmn_grays=np.transpose(np.where((new_subject_labels_refined == 0) | (new_subject_labels_refined == 14) | (new_subject_labels_refined == 15)))
        
        # run template matching to refine dmn_grays
        for i in range(dmn_grays.shape[0]):
            
            if i % 5000 == 0:
                #print('calculating dmn greyordinate: ' + str([i]))
                print('\t\tdefault mode grey refinement has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
            
            # j will loop through each of the candidate networks 
            for j in range(eta_to_template_vox_modified.shape[1]):
                
                # skip networks which don't exist
                if j in (3,5,16):
                    pass
                
                else:
                    
                    # identify voxels which are not 'NaN' in cmap and tmap
                    dmn_gray = dmn_grays[i].item()
                    goodvox = np.where(~np.isnan(corr_mat_full[dmn_gray,:]) & ~np.isnan(cifti_template_mat_full_refined[:,j]))
                    # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                    cmap = corr_mat_full[dmn_gray,(goodvox)]
                    # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                    tmap = cifti_template_mat_full_refined[goodvox,j]
                    
                    # check if cmap or tmap are zero, and if so return 'nan'
                    if (np.sum(cmap) == 0) and (np.sum(tmap) == 0):
                        eta_to_template_vox_modified[dmn_gray,j] = float("nan")
                        r_to_template_vox_modified[dmn_gray,j] = float("nan")
                        
                    # if cmap and tmap are not zero calculate effect size as usual
                    else:
                        
                        # calculation of eta values using RH within and between sum of squares ##########################################
                        
                        # calculate the grand mean across all elements of tmap and cmap
                        Mgrand = (np.mean(tmap) + np.mean(cmap)) / 2
                        # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                        Mwithin = np.mean([tmap, cmap], axis=0)
                        # get a sum of squares within
                        SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                        # get a sum of squares total
                        SStot    = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                        # eta is calculated as 1 - this proportional variance explained
                        eta = (1 - (SSwithin/SStot))
                        eta_to_template_vox_modified[dmn_gray,j] = (1 - (SSwithin/SStot))
                        
                        # calculation of r values using correlation ######################################################################
                        
                        r_to_template_vox_modified[dmn_gray,j] = np.ma.corrcoef(cmap, tmap)[0,1]
                        
                        # delete outputs before the next loop
                        del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox
        
        print('\t\trefinement of default mode grays complete')
        print('\t\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
        print('\t\t')
    
    # refine greyordinates in SAL(8) or CO(9)
    refine_salience_grays=False
    if refine_salience_grays:
        
        print('\trefinement of salience grays:')
        steptimer=time.time()
        
        # find the networks that belong to SAL(8), CO(9), using 7,8 b/c base 0
        sal_grays=np.transpose(np.where((new_subject_labels_refined == 7) | (new_subject_labels_refined == 8)))
        
        # run template matching looping over sal_grays
        for i in range(sal_grays.shape[0]):
            
            if i % 5000 == 0:
                print('\t\tsalience grey refinement has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
            
            # j will loop through each of the candidate networks 
            for j in range(eta_to_template_vox_modified.shape[1]):
                
                # skip networks which don't exist
                if j in (3,5,16):
                    pass
                
                else:
                    
                    # identify voxels which are not 'NaN' in cmap and tmap
                    sal_gray = sal_grays[i].item()
                    goodvox = np.where(~np.isnan(corr_mat_full[sal_gray,:]) & ~np.isnan(cifti_template_mat_full_refined[:,j]))
                    # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                    cmap = corr_mat_full[sal_gray,(goodvox)]
                    # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                    tmap = cifti_template_mat_full_refined[goodvox,j]
                    
                    # check if cmap and tmap are zero, and if so return 'nan'
                    if (np.sum(cmap) == 0) and (np.sum(tmap) == 0):
                        eta_to_template_vox_modified[sal_gray,j] = float("nan")
                        r_to_template_vox_modified[sal_gray,j] = float("nan")
                    
                    # if cmap and tmap are not zero calculate effect size as usual
                    else:
                        
                        # calculation of eta values using RH within and between sum of squares ##########################################
                        
                        # calculate the grand mean across all elements of tmap and cmap
                        Mgrand = (np.mean(tmap) + np.mean(cmap)) / 2
                        # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                        Mwithin = np.mean([tmap, cmap], axis=0)
                        # get a sum of squares within
                        SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                        # get a sum of squares total
                        SStot    = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                        # eta is calculated as 1 - this proportional variance explained
                        # eta = (1 - (SSwithin/SStot))
                        eta_to_template_vox_modified[sal_gray,j] = (1 - (SSwithin/SStot))
                        
                        # calculation of r values using correlation ######################################################################
                        
                        r_to_template_vox_modified[sal_gray,j] = np.ma.corrcoef(cmap, tmap)[0,1]
                        
                        # delete outputs before the next loop
                        del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox
        
        print('\t\trefinement of salience grays complete')
        print('\t\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
    
    print('\t')
    print('\twinner-take-all re-assignment of greyordinates to networks')
    print('\t')
    new_subject_labels_refined = np.nanargmax(eta_to_template_vox_modified,axis=1)
    new_subject_labels_r_refined = np.nanargmax(r_to_template_vox_modified,axis=1)
    
    print('\ttemplate matching refinement network complete')
    print('\tindividual step took: ' + str(round((time.time()-refinetimer)/60,3)) + ' minutes')
    print('\t')

save_outputs=True
if save_outputs:
    
    print('saving outputs')
    print('\t')
    steptimer=time.time()
    
    print('\tsaving intermediates to .mat file')
    
    # create different dictionaries depent on running refine_template_matching
    if refine_template_matching:
        
        # b/c python is base zero, add 1 to each element of new_subject_labels to change network labels from 0:17 -> 1:18
        new_subject_labels = np.array([x + 1 for x in new_subject_labels])
        new_subject_labels_refined = np.array([x + 1 for x in new_subject_labels_refined])
        new_subject_labels_r = np.array([x + 1 for x in new_subject_labels_r])
        new_subject_labels_r_refined = np.array([x + 1 for x in new_subject_labels_r_refined])
        
        # combine outputs into a dictionary
        mat_dict={"network_names": np.array(network_names,dtype=object)
                ,"eta_to_template_vox": eta_to_template_vox
                ,"eta_to_template_vox_modified": eta_to_template_vox_modified
                ,"new_subject_labels": new_subject_labels
                ,"new_subject_labels_refined": new_subject_labels_refined
                ,"r_to_template_vox": r_to_template_vox
                ,"r_to_template_vox_modified": r_to_template_vox_modified
                ,"new_subject_labels_r": new_subject_labels_r
                ,"new_subject_labels_r_refined": new_subject_labels_r_refined
                }
        
        # the actual saving of the dictionary with scipy
        sio.savemat(matfile_out,mat_dict,appendmat=False)
        
    else:
        new_subject_labels = np.array([x + 1 for x in new_subject_labels]).reshape(1,greyordinates)
        new_subject_labels_r = np.array([x + 1 for x in new_subject_labels_r]).reshape(1,greyordinates)
        mat_dict={"network_names": np.array(network_names,dtype=object)
                ,"eta_to_template_vox": eta_to_template_vox
                ,"new_subject_labels":new_subject_labels
                ,"r_to_template_vox":r_to_template_vox
                ,"new_subject_labels_r":new_subject_labels_r}
        # the actual saving of the dictionary with scipy
        sio.savemat(matfile_out,mat_dict,appendmat=False)
    
    print('\tsaving intermediates to .mat file complete')
    print('\t')
    
    print('\tsaving subject dscalar file')
    
    # save the dscalar without refinement network
    dscalar_template = nb.load(dscalar_template_infile)
    scalar_axis = nb.cifti2.ScalarAxis([]) 
    new_img = nb.Cifti2Image(new_subject_labels.reshape(1,greyordinates)
                            ,header=dscalar_template.header
                            ,nifti_header=dscalar_template.nifti_header
                            ,extra=scalar_axis)
    new_img.to_filename(dscalar_outfile)
    
    # save a second dscalar with refinement, if applicable
    if refine_template_matching:
        dscalar_template = nb.load(dscalar_template_infile)
        scalar_axis = nb.cifti2.ScalarAxis([])
        new_img = nb.Cifti2Image(new_subject_labels_refined.reshape(1,greyordinates)
                                ,header=dscalar_template.header
                                ,nifti_header=dscalar_template.nifti_header
                                ,extra=scalar_axis)
        new_img.to_filename(dscalar_outfile_refined)
    
    print('\tsaving subject dscalar file complete')
    print('\t')
    
    # wrap up this step
    print('\tsaving outputs complete')
    print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
    print('\t')