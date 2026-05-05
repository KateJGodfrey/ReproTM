#import glob
import os
import subprocess
import sys
import nibabel as nb
from nibabel import cifti2
import scipy.io as sio       #scipy is for accessing mat files
import numpy as np
import time
import pandas as pd
import mat73                 #mat73 will allow you to load mat files, useful if .mat output already created

# NOTES
# need to test what outputs look like with Z-scored dconn
# should tell users whether or not their input cifti and input template are same units
# need to figure out how to put refined greyordinates back into matrix

#import ciftify as cifti
wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/'
dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn'
cifti_output_folder='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching'
dconn_in='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm_zscoredRH.dconn.nii'
template_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat'
matfile_out='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching2.0._Z/sub-10227_ses-combined_task-cross_TM.mat'
output_cifti_name='sub-10227_ses-combined_task-crossMENORDICtrimmed_smoothed2mm_censoredFD2mm_zscoredRH_TM'
TEMPLATEMINIMUM = 1.00
SCANTEMPLATEMINIMUM = 3.00
wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/'
network_names = ['DMN','Vis','FP','','DAN','','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON','','SCAN']
refine_scan_network=True

# if debug_mode = True, script will print out operations in very verbose style
debug_mode=False

#load the correlation matrix
print('loading the correlation matrix')
dconn_filename = dir_in + dconn_file
subject_cii=nb.load(dconn_filename)
corr_mat_full=subject_cii.get_fdata()   # get correlation matrix out of input cifti
del subject_cii                         # save memory

#load the template and assign seed matrix to variable
print('loading the template: ' + template_path)
print(' ')
cifti_template = sio.loadmat(template_path)
cifti_template_mat_full = cifti_template.get("seed_matrix")

#remove seeds below the threshold
print('template minimum is set at: ' + str(TEMPLATEMINIMUM))
print('seeds below template minimum being set to NaN')
cifti_template_mat_full[cifti_template_mat_full < TEMPLATEMINIMUM] = 'NaN'

# create empty array to store template matched outputs (1 row per greyordinate and 1 column per network)
# eta_to_template_vox = calculate association to template using RH original sum of squares
# r_to_template_vox = calculate association to template with correlation
eta_to_template_vox = np.zeros((int(len(corr_mat_full)),int(len(network_names))))
r_to_template_vox = np.zeros((int(len(corr_mat_full)),int(len(network_names))))

# get the total number of greyordinates to specify number of loops
greyordinates = corr_mat_full.shape[0]
print('number of greyordinates: ' + str(greyordinates))

# start template matching
print('starting calculation of similarity (eta) of verteces to template')
print('template matching started at: ' + time.ctime())
totaltimer = time.time()
run_template_matching = True

if run_template_matching: 

    # i will loop through each vertex from the correlation matrix
    for i in range(greyordinates):

        # print a progress update every 5000 voxels
        if i % 5000 == 0:
            print('calculating greyordinate: ' + str(i))
            print('template matching has run for: ' + str(round(round(time.time()-totaltimer,3)/60,3)) + ' minutes')
        
        # j will loop through each of the candidate networks to match to
        for j in range(len(network_names)):

            # skip networks 4, 6, 17 (which don't exist)
            # specified as 3, 5, 16 b/c python is base 0
            if j in (3,5,16):
                pass

            else:
                if debug_mode:
                    #setting i and j for debugging purposes
                    #setting i = 1 tells it we are looking at first voxel
                    #setting j = 1 tells it we are looking at first network
                    i = 0
                    j = 0

                # identify voxels which are not 'NaN' in cmap and tmap
                # 'NaN' greyordinates were specified earlier using a user-specified tmap threshold
                goodvox = np.where(~np.isnan(corr_mat_full[i,:]) & ~np.isnan(cifti_template_mat_full[:,j]))
                # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                cmap = corr_mat_full[i,goodvox]
                # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                tmap = cifti_template_mat_full[goodvox,j]

                # inspect and write out information about cmap and tmap
                # if debug_mode:

                #     # get info about cmap and tmap
                #     print('cmap shape: ' + str(cmap.shape))
                #     print('first 5 values of cmap: ')
                #     print(cmap[0:5])
                #     print('')
                #     print('tmap shape: ' + str(tmap.shape))
                #     print('first 5 values of tmap: ')
                #     print(tmap[0:5])

                #     # write cmap and tmap to a text file to inspect and use with another software
                #     cmap_outfile = dir_in + '/' + 'sub-10227_ses-combined_task-cross_greyordinate-1_cmap.txt'
                #     np.savetxt(cmap_outfile, cmap, fmt="%.4e")

                #     tmap_outfile = dir_in + '/' + 'sub-10227_ses-combined_task-cross_network-1_tmap.txt'
                #     np.savetxt(tmap_outfile, tmap, fmt="%.4e")

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
                # eta = (1 - (SSwithin/SStot))
                eta_to_template_vox[i,j] = (1 - (SSwithin/SStot))

                # calculation of r values using correlation ######################################################################

                r_to_template_vox[i,j] = np.corrcoef(cmap, tmap)[0,1]

                # delete outputs before the next loop
                del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox

# can also get eta_to_template_vox from a previously created .mat file
mat_filepath='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/sub-10227_ses-combined_task-crossMENORDICtrimmed_TM.mat'
mat_data = mat73.loadmat(mat_filepath)
eta_to_template_vox = pd.DataFrame(mat_data.get('eta_to_template_vox'))

# assign greyordinates to networks with a winner take all approach ###############################

# get new network assignments by using idxmax to return which column has the highest eta
new_subject_labels = np.zeros(int(len(corr_mat_full)))
new_subject_labels = eta_to_template_vox.idxmax(axis=1)

# get new network assignments by using idxmax to return which column has the highest r 
new_subject_labels_r = np.zeros(int(len(corr_mat_full)))
new_subject_labels_r = r_to_template_vox.idxmax(axis=1)

# assign greyordinates to networks with a winner take all approach ###############################

# need to refine the template matching threshold to identify scan network ########################
if refine_scan_network:
    print('template matching for identifying scan network has begun')
    print('scan network threshold: ' + str(SCANTEMPLATEMINIMUM))
    print('')

    # copy winner-take-all-results into an array for scan network
    new_subject_labels_scan=new_subject_labels

    # get the seed matrix from before and but apply the new scan network threshold
    cifti_template_mat_full_scan = cifti_template.get("seed_matrix")
    cifti_template_mat_full_scan[cifti_template_mat_full_scan < SCANTEMPLATEMINIMUM] = 'NaN'

    print('refinement of motor grays has begun')
    print('individual step beginning at: ' + str(round(round(time.time(),3)/60,3)))
    steptimer=time.time()

    # find the networks that belong to SMd(10), SMl(11), or SCAN (18), using 9,10,17 b/c base 0
    motor_grays=np.transpose(np.where((new_subject_labels_scan == 9) | (new_subject_labels_scan == 10) | (new_subject_labels_scan == 17)))

    # create empty array to store template matched outputs at scan threshold for motor grays
    eta_to_template_vox_scan = np.zeros((int(len(motor_grays)),int(len(network_names))))
    r_to_template_vox_scan = np.zeros((int(len(motor_grays)),int(len(network_names))))

    # run template matching looping over all motor_grays and all networks
    for i in range(eta_to_template_vox_scan.shape[0]):

        if i % 5000 == 0:
            print('calculating motor greyordinate: ' + str(motor_grays[i].item()))
            print('template matching for scan network has run for: ' + str(round(round(time.time()-steptimer,3)/60,3)) + ' minutes')
        
        # j will loop through each of the candidate networks 
        for j in range(eta_to_template_vox_scan.shape[1]):

            # skip networks which don't exist
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

                # check if cmap and tmap are zero, and if so return 'nan'
                if (np.sum(cmap) == 0) and (np.sum(tmap) == 0):
                    eta_to_template_vox_scan[i,j] = float("nan")
                    r_to_template_vox_scan[i,j] = float("nan")

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
                    eta_to_template_vox_scan[i,j] = (1 - (SSwithin/SStot))

                    # calculation of r values using correlation ######################################################################

                    r_to_template_vox_scan[i,j] = np.corrcoef(cmap, tmap)[0,1]

                    # delete outputs before the next loop
                    del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox

# NOT TESTED: ###################################################################################
# Kate to test that indexing of grays back into matrix is working properly
# assign new eta values to output

# eta_to_template_vox_modified = eta_to_template_vox
# eta_to_template_vox_modified[motor_grays,] = eta_to_template_vox_scan
# eta_to_template_vox_modified[motor_grays_index==True,] = eta_to_template_vox_scan

# # assign new r values to output
# r_to_template_vox_modified = r_to_template_vox
# r_to_template_vox_modified[motor_grays,] = r_to_template_vox_scan

# NOT TESTED: ###################################################################################

print('refinement of motor grays complete')
print('individual step completed at: ' + str(round(round(time.time(),3)/60,3)))
print('')

print('refinement of default mode grays has begun')
print('individual step beginning at: ' + str(round(round(time.time(),3)/60,3)))
steptimer=time.time()

# find the networks that belong to DMN(1), PMN(15), or PON (16), using 0,14,15 b/c base 0
dmn_grays=np.transpose(np.where((new_subject_labels_scan == 0) | (new_subject_labels_scan == 14) | (new_subject_labels_scan == 15)))

# create empty array to store template matched outputs at scan threshold for default mode grays
eta_to_template_vox_dmn = np.zeros((int(len(dmn_grays)),int(len(network_names))))
r_to_template_vox_dmn = np.zeros((int(len(dmn_grays)),int(len(network_names))))

# run template matching looping over dmn_grays
for i in range(eta_to_template_vox_dmn.shape[0]):

    if i % 5000 == 0:
        print('calculating dmn greyordinate: ' + str(dmn_grays[i].item()))
        print('template matching for dmn grayordinates network has run for: ' + str(round(round(time.time()-steptimer,3)/60,3)) + ' minutes')
    
    # j will loop through each of the candidate networks 
    for j in range(eta_to_template_vox_dmn.shape[1]):

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

            # check if cmap and tmap are zero, and if so return 'nan'
            if (np.sum(cmap) == 0) and (np.sum(tmap) == 0):
                eta_to_template_vox_dmn[i,j] = float("nan")
                r_to_template_vox_dmn[i,j] = float("nan")

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
                eta_to_template_vox_dmn[i,j] = (1 - (SSwithin/SStot))

                # calculation of r values using correlation ######################################################################

                r_to_template_vox_dmn[i,j] = np.corrcoef(cmap, tmap)[0,1]

                # delete outputs before the next loop
                del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox

# NOT TESTED: ###################################################################################
# assign new eta values to output
eta_to_template_vox_modified[dmn_grays,:] = eta_to_template_vox_dmn

# assign new r values to output
r_to_template_vox_modified[dmn_grays,:] = r_to_template_vox_dmn
# NOT TESTED: ###################################################################################

print('refinement of default mode grays complete')
print('individual step completed at: ' + str(round(round(time.time(),3)/60,3)))
print('')

print('refinement of salience grays has begun')
print('individual step beginning at: ' + str(round(round(time.time(),3)/60,3)))
steptimer=time.time()

# find the networks that belong to SAL(8), CO(9), using 7,8 b/c base 0
sal_grays=np.transpose(np.where((new_subject_labels_scan == 7) | (new_subject_labels_scan == 8)))

# create empty array to store template matched outputs at scan threshold for default mode grays
eta_to_template_vox_sal = np.zeros((int(len(sal_grays)),int(len(network_names))))
r_to_template_vox_sal = np.zeros((int(len(sal_grays)),int(len(network_names))))

# run template matching looping over sal_grays
for i in range(eta_to_template_vox_sal.shape[0]):

    if i % 5000 == 0:
        print('calculating sal greyordinate: ' + str(sal_grays[i].item()))
        print('template matching for sal grayordinates network has run for: ' + str(round(round(time.time()-steptimer,3)/60,3)) + ' minutes')
    
    # j will loop through each of the candidate networks 
    for j in range(eta_to_template_vox_sal.shape[1]):

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

            # check if cmap and tmap are zero, and if so return 'nan'
            if (np.sum(cmap) == 0) and (np.sum(tmap) == 0):
                eta_to_template_vox_sal[i,j] = float("nan")
                r_to_template_vox_sal[i,j] = float("nan")

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
                eta_to_template_vox_sal[i,j] = (1 - (SSwithin/SStot))

                # calculation of r values using correlation ######################################################################

                r_to_template_vox_sal[i,j] = np.corrcoef(cmap, tmap)[0,1]

                # delete outputs before the next loop
                del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox

# NOT TESTED: ###################################################################################
# assign new eta values to output
eta_to_template_vox_modified[sal_grays,:] = eta_to_template_vox_sal

# assign new r values to output
r_to_template_vox_modified[sal_grays,:] = r_to_template_vox_sal
# NOT TESTED: ###################################################################################

print('refinement of salience grays complete')
print('individual step completed at: ' + str(round(round(time.time(),3)/60,3)))
print('')

print('refining winner-take-all network assignments')
new_subject_labels_scan = eta_to_template_vox_modified.idxmax(axis=1)

print('saving template matching outputs to .mat file')
mat_dict={"eta_to_template_vox": eta_to_template_vox
        ,"r_to_template_vox":r_to_template_vox
        ,"new_subject_labels":new_subject_labels}

if refine_scan_network:  #overwrite to include 
    mat_dict={"eta_to_template_vox": eta_to_template_vox
            ,"r_to_template_vox":r_to_template_vox
            ,"new_subject_labels":new_subject_labels
            ,"new_subject_labels_scan":new_subject_labels_scan}


## WILL NEED TO ADD ONE TO THE NETWORK ASSIGNMENTS AT END OF SCRIPT B/C PYTHON IS BASE 0

   