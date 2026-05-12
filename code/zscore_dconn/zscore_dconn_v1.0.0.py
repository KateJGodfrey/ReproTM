# this script converts CIFTI dconn.nii to a CIFTI Z-scored dconn.nii

# the input whole-brain dconn correlation matrix will be Z-scored separately for:
# i.   each hemisphere
# ii.  the subcortical region
# iii. the connections between the cortex and the subcortex 
# sectioning the dconn normalizes connectivity between subcortex and cortex when considering potential for decreased SNR in subcortex
# source code and Z-scoring method described above taken from Hermisillo et al (2024) (https://doi.org/10.1038/s41593-024-01596-5)

# author: Kate J. Godfrey 
# github: https://github.com/KateJGodfrey/ReproTM

# load packages
import numpy as np
import nibabel as nb
from scipy.stats import zscore
import argparse
import textwrap

parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter
                             ,description=textwrap.dedent('''Transform dconn to a Z-scored dconn\ntransformation performed separately by cortical hemisphere and cortex/subcortex.'''))
parser.add_argument('--dconn_infile',required=True,help="path to input dconn")
parser.add_argument('--dconn_outfile',required=True,help="path to output z-scored dconn")
args = parser.parse_args()

dconn_infile=args.dconn_infile
dconn_outfile=args.dconn_outfile

print('')
print('input dconn:')
print(dconn_infile)
print('')

print('output z-scored dconn:')
print(dconn_outfile)
print('')

# uncomment to hardcode necessary paths, useful for debugging
# wb_command='/path/workbench'
# dconn_infile='/path/input_dconn.dconn.nii'
# dconn_outfile='/path/output_dconn_stat-zscored.dconn.nii'

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
