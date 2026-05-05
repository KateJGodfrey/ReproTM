# this script will take either a Fisher Z or Pearson correlation dconn and transform to a Z score matrix
# the input dconn correlation matrix of cortical and subcortical regions will be Z-scored separately for:
# i.   each hemisphere
# ii.  the subcortical region
# iii. the connections between the cortex and the subcortex 
# sectioning the dconn normalizes connectivity between subcortex and cortex when considering a potential of decreased SNR in subcortex
# this Z-scoring method was originally described in Hermisillo et al (2024) (https://doi.org/10.1038/s41593-024-01596-5)

# source matlab code for this script can be found at:
# /home/faird/shared/code/internal/utilities/Zscore_dconn/Zscore_dconn.m

# load packages
import numpy as np
import nibabel as nb
from scipy.stats import zscore
import sys

# set paths
wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/'
dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn'
dconn_infile='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm.dconn.nii'
dconn_outfile='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm_zscoredRH.dconn.nii'

# load data
print('loading the correlation matrix')
dconn_infile_path = dir_in + dconn_infile
subject_cii=nb.load(dconn_infile_path)
dconn=subject_cii.get_fdata()   # get correlation matrix out of input cifti

debug=False
if debug:
    # can do some calculations on dconn
    num_nans = np.count_nonzero(np.isnan(dconn))   # number of NaN values
    matrix_sum = np.nansum(dconn)                  # sum of all elements that aren't NaN
    print('dconn size: ' + str(dconn.size))
    print('dconn number of NaN: ' + str(num_nans))
    print('dconn sum of all elements: ' + str(matrix_sum))
    del num_nans, matrix_sum

# section dconn by hemisphere and cortical/ subcortical connections
# note: difference between python and matlab script indexing:
# python starts at base 0, while matlab is base 1
# python indexing goes up to, but doesn't include, the second argument
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

del dconn # save memory

if debug:
    # can compare matrix calculated here to the matlab code
    test_matrix = SS
    num_nans = np.count_nonzero(np.isnan(test_matrix))   # number of NaN values
    matrix_sum = np.nansum(test_matrix)
    test_matrix.size
    print('number of NaN: ' + str(num_nans))
    print('matrix sum: ' + str(matrix_sum))
    print('matrix size: ' + str(test_matrix.size))

# Z-score across all elements of sectioned matrices 
# NaN values are ignored and returns 'NaN' to those elements
print('zscore beginning')
ZLL_mat = zscore(LL,axis=None,nan_policy='omit')  # axis=None: z-score includes all elements for mean & std
ZLR_mat = zscore(LR,axis=None,nan_policy='omit')  # nan_policy='omit': ignore nan when calculating z-score
ZLS_mat = zscore(LS,axis=None,nan_policy='omit')  #                    and return nan back to nan elements 
ZRL_mat = zscore(RL,axis=None,nan_policy='omit')
ZRR_mat = zscore(RR,axis=None,nan_policy='omit')
ZRS_mat = zscore(RS,axis=None,nan_policy='omit')
ZSL_mat = zscore(SL,axis=None,nan_policy='omit')
ZSR_mat = zscore(SR,axis=None,nan_policy='omit')
ZSS_mat = zscore(SS,axis=None,nan_policy='omit')
print('zscore complete')

# can compare matrix calculated here to the matlab code
test_matrix = ZLL_mat
num_nans = np.count_nonzero(np.isnan(test_matrix))   # number of NaN values
matrix_sum = np.nansum(test_matrix)
test_matrix.size
print('test_matrix: ZLLmat')
print('number of NaN: ' + str(num_nans))
print('matrix sum: ' + str(matrix_sum))
print('matrix size: ' + str(test_matrix.size))

# rewrite matrix
print('rewriting matrix beginning')
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
print('rewriting matrix complete')

del ZLL_mat, ZLR_mat, ZLS_mat, ZRL_mat, ZRR_mat, ZRS_mat, ZSL_mat, ZSR_mat, ZSS_mat

# give the new dconn values to the subjects cii and save with new name
print('writing output dconn beginning')
subject_cii.cdata = newdconn
dconn_outfile_path = dir_in + dconn_outfile
nb.save(subject_cii, dconn_outfile_path)

if debug:
    x = np.array([[25.11, 30.10, np.nan, 32.02, 43.15],
                [14.95, 16.06, 121.25, 94.35, 29.81]])
    y = [25.11, 30.10, 32.02, 43.15, 14.95, 16.06, 121.25, 94.35, 29.81]
    stats.zscore(x, axis=None,nan_policy='omit')