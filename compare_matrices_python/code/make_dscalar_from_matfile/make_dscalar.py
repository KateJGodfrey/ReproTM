import os
import subprocess
import sys
import nibabel as nb
from nibabel import cifti2
import scipy.io as sio       # scipy is for accessing mat files
import numpy as np
import time
import pandas as pd

# make a cifti dscalar network map from a matfile

# paths can also be hardcoded to help with debugging
# wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/'
# cifti_output_folder = '/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching'
# matfile_in = '/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_fisherZ_zscoreRH.dconn.nii'
# dscalar_outfile ='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/v6_interpolated/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_fisherZ_zscoreRH_pythonTM.dscalar.nii'
# dir_in = '/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching'
# mat_filepath='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/sub-11421_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.mat'
# dscalar_outfile ='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/sub-11421_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonScanTM.dscalar.nii'
# dscalar_template_infile ='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/support_files/sub-10227_ses-combined_task-cross_grims_recoloredCifti2.dscalar.nii'

# input template and template options
mat_filepath = [sys.argv[1]][0]
dscalar_outfile=[sys.argv[2]][0]
dscalar_template_infile = [sys.argv[3]][0]

print('loading mat file')
mat_data = sio.loadmat(mat_filepath)
eta_to_template_vox = np.array(mat_data.get('eta_to_template_vox'))
eta_to_template_vox_modified = np.array(mat_data.get('eta_to_template_vox_scan'))
new_subject_labels = np.array(mat_data.get('new_subject_labels'))
new_subject_labels_scan = np.array(mat_data.get('new_subject_labels_scan'))

print('making dscalar')
dscalar_template = nb.load(dscalar_template_infile)
header = dscalar_template.header
spatial_axes = header.get_axis(1)
scalar_axis = nb.cifti2.ScalarAxis([])

new_img = nb.Cifti2Image(new_subject_labels_scan.reshape(1,new_subject_labels_scan.shape[1])
                        ,header=dscalar_template.header
                        ,nifti_header=dscalar_template.nifti_header
                        ,extra=scalar_axis)
new_img.to_filename(dscalar_outfile)

# wrap up this step
print('\t')
print('\tsaving outputs complete')
print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
print('\t')