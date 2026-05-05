#!/bin/bash -l

#SBATCH -J dconnZ
#SBATCH -c 1
#SBATCH --tmp=200gb
#SBATCH --mem=200gb
#SBATCH -t 1:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kgodfrey@umn.edu
#SBATCH -p msibigmem
#SBATCH -o output_logs/dconnZ_%j.out
#SBATCH -e output_logs/dconnZ_%j.err
#SBATCH -A btervocl

# user specified terminal inputs
SUB=${1}
SES=${2}
TASK=${3}

# update output log with script details
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "job: ${SLURM_JOB_ID}"
echo "sub: sub-${SUB}"
echo "ses: ses-${SES}"
echo "task: task-${TASK}"
echo " "

# set input directory
dir_in=/home/btervocl/shared/projects/MINT/kgodfrey/template_matching

# set paths for input and output dconns
dconn_in=${dir_in}/derivatives/dconn/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm.dconn.nii
dconn_out=${dir_in}/derivatives/dconn/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH.dconn.nii
dconn_out=/scratch.global/kgodfrey/zscore_dconn/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH.dconn.nii

# link default output log to new log with subject identifier
ln -f ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.out ${dir_in}/output_logs/dconnZ_sub-${SUB}_ses-${SES}_task-${TASK}.out
ln -f ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.err ${dir_in}/output_logs/dconnZ_sub-${SUB}_ses-${SES}_task-${TASK}.err

# source python and run script to z-score dconns
source /home/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
code_path=/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/zscore_dconn/Zscore_dconn_v2.1.py
python -u ${code_path} ${dconn_in} ${dconn_out}

echo ""
echo "script completed at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"

# remove original output logs
rm ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.out
rm ${dir_in}/output_logs/dconnZ_${SLURM_JOB_ID}.err