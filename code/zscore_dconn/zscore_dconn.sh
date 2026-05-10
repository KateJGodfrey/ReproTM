#!/bin/bash -l

#SBATCH -J dconnZ
#SBATCH -c 1
#SBATCH --tmp=200gb
#SBATCH --mem=200gb
#SBATCH -t 1:00:00
#SBATCH -p your_partition
#SBATCH -o dconnZ_%j.out
#SBATCH -e dconnZ_%j.err

# author: Kate J. Godfrey 
# email:  godfreykatej@gmail.com 
# github: https://github.com/KateJGodfrey/ReproTM

# this script converts CIFTI dconn.nii to a CIFTI Z-scored dconn.nii

# the input dconn of cortical and subcortical regions will be Z-scored separately for:
# i.   each hemisphere
# ii.  the subcortical region
# iii. the connections between the cortex and the subcortex 
# dconn sectioning normalizes connectivity between subcortex and cortex when considering potential for decreased SNR in subcortex
# source code and Z-scoring method described above taken from Hermisillo et al (2024) (https://doi.org/10.1038/s41593-024-01596-5)

# update output log with script details
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"

# user specified terminal inputs
SUB=${1}
SES=${2}
TASK=${3}

# set directories and paths
dir_in=/path/to/your/data/directory
dconn_indir=${dir_in}/derivatives
dconn_outdir=${dir_in}/derivatives

# set paths for input and output dconns
dconn_infile=${dconn_indir}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}.dconn.nii
dconn_outfile=${dconn_outdir}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_stat-zscored.dconn.nii

# set path for python code 
code_path=/*/ReproTM/zscore_dconn/zscore_dconn_v1.0.0.py

# check if input dconn exists
if [ -e ${dconn_in} ]
then
    echo "input dconn file exists"
    echo "running script to z-score"
    echo " " 
    
    # make sure output folder exists, if not create it
    if [ ! -d ${dconn_outdir}/sub-${SUB}/ses-${SES}/func ]; then
        mkdir -p ${dconn_outdir}/sub-${SUB}/ses-${SES}/func
    fi  

    # run script to z-score dconns
    python -u ${code_path} --dconn_infile ${dconn_infile} --dconn_outfile ${dconn_outfile}
    
else
    echo "input dconn file does not exist"
    echo "not running script to z-score"
    echo " " 
fi

echo ""
echo "script completed at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"