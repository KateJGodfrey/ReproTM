# this script will take either a Fisher Z or Pearson correlation dconn and transform to a Z score dconn
# the input dconn correlation matrix of cortical and subcortical regions will be Z-scored separately for:
# i.   each hemisphere
# ii.  the subcortical region
# iii. the connections between the cortex and the subcortex 
# sectioning the dconn normalizes connectivity between subcortex and cortex when considering potential for decreased SNR in subcortex
# the Z-scoring method described above taken from Hermisillo et al (2024) (https://doi.org/10.1038/s41593-024-01596-5)

# source matlab code for the python script in */src or at path:
# /home/faird/shared/code/internal/utilities/Zscore_dconn/Zscore_dconn.m

# update log
# v1.0 (major update)
# i. updated hardcode script to run with terminal command line arguements
# v2.0 (major update)
# i.  update to properly save the z-scored dconn (previously was writing back the old dconn)
# v2.1 (minor update)
# i.  updated what information is printed with 'debug' option
# ii. cleaned up old code for readability

# load packages
import numpy as np
import nibabel as nb
from scipy.stats import zscore
import sys

dconn_infile = [sys.argv[1]][0]
dconn_outfile = [sys.argv[2]][0]

print('input dconn:')
print(dconn_infile)
print('')

print('output z-scored dconn:')
print(dconn_outfile)
print('')

# set paths (if you want to run hardcoded script from command line)
# wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/'
# dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn'
# dconn_infile='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm.dconn.nii'
# dconn_outfile='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm_zscoredRH.dconn.nii'

print('loading the correlation matrix')
subject_cii=nb.load(dconn_infile)
dconn=subject_cii.get_fdata()   # get correlation matrix out of input cifti

debug=True
if debug:
    print('first few elements of dconn: ')
    print(dconn[0:4,0:4])

# section dconn by hemisphere and cortical/ subcortical connections
# note: differences between python and matlab script indexing:
# i.  python starts at base 0, while matlab is base 1
# ii. python indexing goes up to (but doesn't include) the second argument
print('dconn sectioning beginning')
LL = dconn[0:29696,0:29696]
LR = dconn[29696:59412,0:29696]
LS = dconn[59412:,0:29696]
RL = dconn[0:29696,29696:59412]
RR = dconn[29696:59412,29696:59412]
RS = dconn[59412:,29696:59412]
SL = dconn[0:29696,59412:]
SR = dconn[29696:59412,59412:]
SS = dconn[59412:,59412:]
print('dconn sectioning complete')

# z-score across all elements of sectioned dconn
# nan values are ignored and 'NaN' is returned back to those elements in Z-score matrix
print('zscore beginning')
ZLL_mat = zscore(LL,axis=None,nan_policy='omit')  # axis=None: z-score includes all elements to calculate mean & std
ZLR_mat = zscore(LR,axis=None,nan_policy='omit')  # nan_policy='omit': ignore nan when calculating z-score
ZLS_mat = zscore(LS,axis=None,nan_policy='omit')  #                    and return nan back to nan elements 
ZRL_mat = zscore(RL,axis=None,nan_policy='omit')
ZRR_mat = zscore(RR,axis=None,nan_policy='omit')
ZRS_mat = zscore(RS,axis=None,nan_policy='omit')
ZSL_mat = zscore(SL,axis=None,nan_policy='omit')
ZSR_mat = zscore(SR,axis=None,nan_policy='omit')
ZSS_mat = zscore(SS,axis=None,nan_policy='omit')
print('zscore complete')

# rewrite matrix
print('re-writing matrix beginning')
newdconn = np.zeros((91282,91282))
newdconn[0:29696,0:29696] = ZLL_mat
newdconn[29696:59412,0:29696] = ZLR_mat
newdconn[59412:,0:29696] = ZLS_mat
newdconn[0:29696,29696:59412] = ZRL_mat
newdconn[29696:59412,29696:59412] = ZRR_mat
newdconn[59412:,29696:59412] = ZRS_mat
newdconn[0:29696,59412:] = ZSL_mat
newdconn[29696:59412,59412:] = ZSR_mat
newdconn[59412:,59412:] = ZSS_mat
print('re-writing matrix complete')

if debug:
    print('first few elements of z-scored dconn: ')
    print(newdconn[0:4,0:4])

# give the z-scored dconn values to subject cii and save with new name cifti name
print('saving output file beginning')
cifti_template = nb.load(dconn_infile)
new_img = nb.Cifti2Image(newdconn, header=cifti_template.header,
                         nifti_header=cifti_template.nifti_header)
new_img.to_filename(dconn_outfile)
print('saving output file complete')