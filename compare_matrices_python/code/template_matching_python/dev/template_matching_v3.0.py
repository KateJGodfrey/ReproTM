import os
import subprocess
import sys
import nibabel as nb
from nibabel import cifti2
import scipy.io as sio       # scipy is for accessing mat files
import numpy as np
import time
import pandas as pd
import mat73                 # mat73 will allow you to load mat files, useful if .mat output already created

# NOTES
# need to test what outputs look like with Z-scored dconn
# should tell users whether or not their input cifti and input template are same units
# does this work without wb_command?

# update log
# v3.0 
# major updates:
# i.   updated eta_to_vox_modified = eta_to_template_vox.copy() to prevent original from overwrite
# ii.  rather than returning NaN if cmap == 0 or tmap == 0 now only returns NaN if SStot = 0
# iii. tried new method of writing dscalar file

# get the current time
totaltimer = time.time()

# set paths
wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/'
dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn'
cifti_output_folder='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching'
# input files
dconn_infile='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm_zscoredRH.dconn.nii'
template_infile='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat'
dscalar_template_infile='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/91282_GreyordinatesCifti2.dscalar.nii'

# output files
matfile_out='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matchingv2.0_Z/sub-10227_ses-combined_task-cross_TM.mat'
dscalar_outfile='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matchingv2.0_Z/sub-10227_ses-combined_task-crossMENORDICtrimmed_smoothed2mm_censoredFD2mm_zscoredRH_TM.dscalar.nii'

# user specified options
TEMPLATEMINIMUM = 1.00
SCANTEMPLATEMINIMUM = 3.00
network_names = ['DMN','Vis','FP','','DAN','','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON','','SCAN']
refine_scan_network=True
debug_mode=False # if debug_mode = True, script will print out operations in very verbose style


file_checker=False
if file_checker:

    # check for input dconn
    if os.path.isfile(dconn_infile) == False:
        print('your dconn path doesnt exist')
        print('\tinput dconn file path: ' + dconn_infile)
    # check for input template
    if os.path.isfile(template_infile) == False:
        print('your template path doesnt exist')
        print('\tinput template file path: ' + template_infile)

load_data=True
if load_data:
    print('loading data')
    print('\t')
    steptimer = time.time()

    print('\tloading the correlation matrix')
    subject_cii=nb.load(dir_in + dconn_infile)
    corr_mat_full=subject_cii.get_fdata()  # get dconn out of input cifti
    greyordinates = corr_mat_full.shape[0]
    del subject_cii # save memory

    print('\tnumber of greyordinates: ' + str(greyordinates))
    print('\t ')

    print('\tloading the template')
    print('\ttemplate minimum set at: ' + str(TEMPLATEMINIMUM))
    print('\ttemplate greyordinates below minimum set to NaN')
    print('\t')
    cifti_template = sio.loadmat(template_infile)
    cifti_template_mat_full = cifti_template.get("seed_matrix")
    cifti_template_mat_full[cifti_template_mat_full < TEMPLATEMINIMUM] = 'NaN'

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
                goodvox = np.where(~np.isnan(corr_mat_full[i,:]) & ~np.isnan(cifti_template_mat_full[:,j]))
                # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                cmap = corr_mat_full[i,goodvox]
                # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                tmap = cifti_template_mat_full[goodvox,j]
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

save_first_step_mat=False
if save_first_step_mat:
    mat_dict={"network_names": network_names
            ,"eta_to_template_vox": eta_to_template_vox
            ,"new_subject_labels":new_subject_labels}
    matfile_out='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matchingv2.0_Z/sub-10227_ses-combined_task-cross_mat-before-refinement.mat'    
    sio.savemat(matfile_out,mat_dict,appendmat=False)

load_old_mat = False
if load_old_mat:
    print('loading old mat file')
    mat_filepath='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matchingv2.0_Z/outputs_v1/sub-10227_ses-combined_task-cross_mat-before-refinement.mat'
    mat_data = sio.loadmat(mat_filepath)
    eta_to_template_vox = np.array(mat_data.get('eta_to_template_vox'))
    new_subject_labels = np.array(mat_data.get('new_subject_labels'))

# need to refine the template matching threshold to identify scan network
refine_scan_network=True 
if refine_scan_network:

    print('template matching for scan network at threshold: ' + str(SCANTEMPLATEMINIMUM))
    print('')
    scantimer = time.time()

    # copy eta results into a new array
    eta_to_template_vox_modified = eta_to_template_vox.copy()
    r_to_template_vox_modified = r_to_template_vox.copy()

    # copy winner-take-all results into new array
    new_subject_labels_scan=new_subject_labels.copy()

    # get the seed matrix from before and but apply the new scan network threshold
    cifti_template_mat_full_scan = cifti_template.get("seed_matrix")
    cifti_template_mat_full_scan[cifti_template_mat_full_scan < SCANTEMPLATEMINIMUM] = 'NaN'

    # refine greyordinates in SMd(10), SMl(11), or SCAN (18) networks
    refine_motor_grays=True
    if refine_motor_grays:
        print('\trefinement of motor grays:')
        #print('individual step beginning at: ' + time.ctime())
        steptimer=time.time()
        # find the greyordinates which belong to SMd(10), SMl(11), or SCAN (18), b/c python is base 0 index w/ 9,10,17
        motor_grays=np.transpose(np.where((new_subject_labels_scan == 9) | (new_subject_labels_scan == 10) | (new_subject_labels_scan == 17)))
        # run template matching looping over all motor_grays and all networks
        for i in range(motor_grays.shape[0]):
            if i % 5000 == 0:
                print('\t\tmotor grey refinement has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
            for j in range(eta_to_template_vox_modified.shape[1]):
                if j in (3,5,16):
                    pass
                else:
                    # identify voxels which are not 'NaN' in cmap and tmap
                    motor_gray = motor_grays[i].item()
                    goodvox = np.where(~np.isnan(corr_mat_full[motor_gray,:]) & ~np.isnan(cifti_template_mat_full_scan[:,j]))
                    # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                    cmap = corr_mat_full[motor_gray,(goodvox)]
                    # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                    tmap = cifti_template_mat_full_scan[goodvox,j]

                    # calculation of eta values using RH within and between sum of squares ##########################################
                    # calculate the grand mean across all elements of tmap and cmap
                    Mgrand = (np.mean(tmap) + np.mean(cmap)) / 2
                    # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                    Mwithin = np.mean([tmap, cmap], axis=0)
                    # get a sum of squares within
                    SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                    # get a sum of squares total
                    SStot    = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                    if(SStot == 0):
                        eta_to_template_vox_modified[motor_gray,j] == float("nan")
                    else:
                        # eta is calculated as 1 - this proportional variance explained
                        eta = (1 - (SSwithin/SStot))
                        eta_to_template_vox_modified[motor_gray,j] = (1 - (SSwithin/SStot))
                    # calculation of r values using correlation ######################################################################
                    r_to_template_vox_modified[motor_gray,j] = np.corrcoef(cmap, tmap)[0,1]
                    # delete outputs before the next loop
                    del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox

        print('\t\trefinement of motor grays complete')
        print('\t\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
        print('\t\t')

    # refine greyordinates in DMN(1), PMN(15), or PON(16)
    refine_default_grays=True
    if refine_default_grays:
        print('\trefinement of default mode grays:')
        steptimer=time.time()

        # find the networks that belong to DMN(1), PMN(15), or PON (16), using 0,14,15 b/c base 0
        dmn_grays=np.transpose(np.where((new_subject_labels_scan == 0) | (new_subject_labels_scan == 14) | (new_subject_labels_scan == 15)))

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
                    goodvox = np.where(~np.isnan(corr_mat_full[dmn_gray,:]) & ~np.isnan(cifti_template_mat_full_scan[:,j]))
                    # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                    cmap = corr_mat_full[dmn_gray,(goodvox)]
                    # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                    tmap = cifti_template_mat_full_scan[goodvox,j]
                    # calculation of eta values using RH within and between sum of squares ##########################################

                    # calculate the grand mean across all elements of tmap and cmap
                    Mgrand = (np.mean(tmap) + np.mean(cmap)) / 2
                    # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                    Mwithin = np.mean([tmap, cmap], axis=0)
                    # get a sum of squares within
                    SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                    # get a sum of squares total
                    SStot = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                    if(SStot == 0):
                        eta_to_template_vox_modified[dmn_gray,j] == float("nan")
                    else:
                        # eta is calculated as 1 - this proportional variance explained
                        #eta = (1 - (SSwithin/SStot))
                        eta_to_template_vox_modified[motor_gray,j] = (1 - (SSwithin/SStot))

                    # calculation of r values using correlation ######################################################################

                    r_to_template_vox_modified[dmn_gray,j] = np.corrcoef(cmap, tmap)[0,1]

                    # delete outputs before the next loop
                    del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox

        print('\t\trefinement of default mode grays complete')
        print('\t\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
        print('\t\t')

    # refine greyordinates in SAL(8) or CO(9)
    refine_salience_grays=True
    if refine_salience_grays:
        print('\trefinement of salience grays:')
        steptimer=time.time()

        # find the networks that belong to SAL(8), CO(9), using 7,8 b/c base 0
        sal_grays=np.transpose(np.where((new_subject_labels_scan == 7) | (new_subject_labels_scan == 8)))

        # run template matching looping over sal_grays
        for i in range(sal_grays.shape[0]):

            if i % 5000 == 0:
                #print('calculating salience greyordinate: ' + str(i))
                print('\t\tsalience grey refinement has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
            
            # j will loop through each of the candidate networks 
            for j in range(eta_to_template_vox_modified.shape[1]):

                # skip networks which don't exist
                if j in (3,5,16):
                    pass

                else:

                    # identify voxels which are not 'NaN' in cmap and tmap
                    sal_gray = sal_grays[i].item()
                    goodvox = np.where(~np.isnan(corr_mat_full[sal_gray,:]) & ~np.isnan(cifti_template_mat_full_scan[:,j]))
                    # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                    cmap = corr_mat_full[sal_gray,(goodvox)]
                    # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                    tmap = cifti_template_mat_full_scan[goodvox,j]
                    # calculation of eta values using RH within and between sum of squares ##########################################

                    # calculate the grand mean across all elements of tmap and cmap
                    Mgrand = (np.mean(tmap) + np.mean(cmap)) / 2
                    # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                    Mwithin = np.mean([tmap, cmap], axis=0)
                    # get a sum of squares within
                    SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                    # get a sum of squares total
                    SStot = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                    # prevent code failure if divide by zero
                    if(SStot == 0):
                        eta_to_template_vox_modified[sal_gray,j] == float("nan")
                    else:
                        # eta is calculated as 1 - this proportional variance explained
                        eta = (1 - (SSwithin/SStot))
                        eta_to_template_vox_modified[sal_gray,j] = (1 - (SSwithin/SStot))

                    # calculation of r values using correlation ######################################################################

                    r_to_template_vox_modified[sal_gray,j] = np.corrcoef(cmap, tmap)[0,1]

                    # delete outputs before the next loop
                    del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox

        print('\t\trefinement of salience grays complete')
        print('\t\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')

    print('\t')
    print('\twinner-take-all re-assignment of greyordinates to networks')
    print('\t')
    new_subject_labels_scan = np.nanargmax(eta_to_template_vox_modified,axis=1)
    new_subject_labels_r_scan = np.nanargmax(r_to_template_vox_modified,axis=1)

    print('\ttemplate matching for scan network refinement complete')
    print('\tindividual step took: ' + str(round((time.time()-scantimer)/60,3)) + ' minutes')
    print('\t')
    
# # can also get eta_to_template_vox from a previously created .mat file
load_old_mat = False
if load_old_mat:
    print('loading old mat file')
    mat_filepath='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matchingv2.0_Z/sub-10227_ses-combined_task-cross_TM.mat'
    mat_data = sio.loadmat(mat_filepath)
    eta_to_template_vox = np.array(mat_data.get('eta_to_template_vox'))
    eta_to_template_vox_modified = np.array(mat_data.get('eta_to_template_vox_scan'))
    new_subject_labels = np.array(mat_data.get('new_subject_labels'))
    new_subject_labels_scan = np.array(mat_data.get('new_subject_labels_scan'))

save_outputs=True
if save_outputs:
    print('saving outputs')
    print('\t')
    steptimer=time.time()

    print('\tsaving intermediates to mat file')
    print('\t')

    # combine outputs into a dictionary
    if refine_scan_network:

        # b/c python is base zero, add 1 to each element of new_subject_labels to number networks from 1:18
        new_subject_labels = np.transpose(np.array([x + 1 for x in new_subject_labels]))
        new_subject_labels_scan = np.transpose(np.array([x + 1 for x in new_subject_labels_scan]))
        new_subject_labels_r = np.transpose(np.array([x + 1 for x in new_subject_labels_r]))
        new_subject_labels_r_scan = np.transpose(np.array([x + 1 for x in new_subject_labels_r_scan]))

        mat_dict={"network_names": network_names
                ,"eta_to_template_vox": eta_to_template_vox
                ,"eta_to_template_vox_modified": eta_to_template_vox_modified
                ,"new_subject_labels": new_subject_labels
                ,"new_subject_labels_scan": new_subject_labels_scan
                ,"r_to_template_vox": r_to_template_vox
                ,"r_to_template_vox_scan": r_to_template_vox_modified
                ,"new_subject_labels_r": new_subject_labels_r
                ,"new_subject_labels_r_scan": new_subject_labels_r_scan
                }
    else:
        new_subject_labels = np.array([x + 1 for x in new_subject_labels])
        new_subject_labels_r = np.array([x + 1 for x in new_subject_labels_r])
        mat_dict={"network_names": network_names
                ,"eta_to_template_vox": eta_to_template_vox
                ,"new_subject_labels":new_subject_labels
                ,"r_to_template_vox":r_to_template_vox
                ,"new_subject_labels_r":new_subject_labels_r}
    
    # save dictionary with scipy
    sio.savemat(matfile_out,mat_dict,appendmat=False)

    print('\tsaving network assignments to dscalar file')
    
    # use a previous dscalar as a template
    dscalar_template=nb.load(dscalar_template_infile)

    # assign new data
    # dscalar_template.cdata = new_subject_labels_scan

    # Create a new Cifti2Image with the new data
    new_cifti = cifti2.Cifti2Image.from_array(new_subject_labels_scan, dscalar_template.header)

    dscalar_template = nb.load(dscalar_template_infile)
    new_img = nb.Cifti2Image(new_subject_labels_scan, header=dscalar_template.header,
                            nifti_header=dscalar_template.nifti_header)
    new_img.to_filename(dscalar_outfile)

    # wrap up this step
    print('\t')
    print('\tsaving outputs complete')
    print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
    print('\t')