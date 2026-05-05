#!/bin/bash -l

#SBATCH -J dconnZ
#SBATCH -c 1
#SBATCH --tmp=200gb
#SBATCH --mem=200gb
#SBATCH -t 1:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pandh015@umn.edu
#SBATCH -p msibigmem
#SBATCH -o /projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/output_logs/dconnZ_%j.out
#SBATCH -e /projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/output_logs/dconnZ_%j.err
#SBATCH -A midb_abcd

# this script runs python code to take either a Fisher Z or Pearson correlation dconn and transform to a Z score dconn
# the input dconn correlation matrix of cortical and subcortical regions will be Z-scored separately for:
# i.   each hemisphere
# ii.  the subcortical region
# iii. the connections between the cortex and the subcortex 
# sectioning the dconn normalizes connectivity between subcortex and cortex when considering potential for decreased SNR in subcortex
# the Z-scoring method described above taken from Hermisillo et al (2024) (https://doi.org/10.1038/s41593-024-01596-5)

# source matlab code for the python script in */code/Zscore_dconn/src or at path:
# /home/faird/shared/code/internal/utilities/Zscore_dconn/Zscore_dconn.m

# update log
# v1.1 (minor updates)
# i.   check that input dconns exist before running
# ii.  user can specify different dconn input and output folders
# iii. script will check for or create a subject-specific folder in output directory

# user specified terminal inputs
SUB=0A4P0LWM
SES=00A
TASK=rest

# update output log with script details
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "job: ${SLURM_JOB_ID}"
echo "sub: sub-${SUB}"
echo "ses: ses-${SES}"
echo "task: task-${TASK}"
echo " "

# set directories and paths
dir_in=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python
dconn_indir=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/data/derivatives/xcp_d_v0.12.0_unstable/cifti_connectivity_outputs
dconn_outdir=/scratch.global/pandh015/precision_maps_via_template_matching/dconnZ
code_path=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/Zscore_dconn/Zscore_dconn_v2.1.py

# set paths for input and output dconns
dconn_in=${dconn_indir}/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_task-${TASK}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii_10_minutes_of_data_at_FD_0.2.dconn.nii
dconn_out=${dconn_outdir}/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_task-${TASK}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii_10_minutes_of_data_at_FD_0.2_Zscore.dconn.nii

# link default output log to new log with subject identifier
ln -f ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.out ${dir_in}/output_logs/dconnZ_sub-${SUB}_ses-${SES}_task-${TASK}.out
ln -f ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.err ${dir_in}/output_logs/dconnZ_sub-${SUB}_ses-${SES}_task-${TASK}.err

# check if input dconn exists
if [ -e ${dconn_in} ]
then
    echo "your input dconn file exists!"
    echo "running script to z-score"
    echo " " 
    
    # first make sure output folder exists, if not create it
    if [ ! -d ${dconn_outdir}/sub-${SUB}/ses-${SES}/func ]; then
        mkdir -p ${dconn_outdir}/sub-${SUB}/ses-${SES}/func
    fi  

    # run script to z-score dconns
    source /projects/standard/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
    python -u ${code_path} ${dconn_in} ${dconn_out}
    
else
    echo "your input dconn file does not exist"
    echo "not running script to z-score"
    echo " " 
fi

echo ""
echo "script completed at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"

# remove original output logs
rm ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.out
rm ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.err