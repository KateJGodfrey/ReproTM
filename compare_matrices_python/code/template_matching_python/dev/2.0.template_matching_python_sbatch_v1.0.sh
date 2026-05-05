#!/bin/bash -l

#SBATCH -J pythonTM
#SBATCH --ntasks=2
#SBATCH --tmp=120gb
#SBATCH --mem=120gb
#SBATCH -t 2:00:00
#SBATCH -p msibigmem
#SBATCH -o /home/faird/shared/code/internal/analytics/compare_matrices_python/output_logs/pythonTM_%j.out
#SBATCH -e /home/faird/shared/code/internal/analytics/compare_matrices_python/output_logs/pythonTM_%j.err
#SBATCH -A btervocl

SUB=${1}
SES=${2}
TASK=${3}

# example command:
# sbatch ./2.0.template_matching_python_sbatch.sh 10227 combined cross

# what directory are you in?
dir_in=/home/faird/shared/code/internal/analytics/compare_matrices_python

# what directory are your input dconns in?
input_folder='/scratch.global/kgodfrey/DCAN_WIP/dconnZ'
# where do you want to save the template matching outputs?
output_folder=${dir_in}/derivatives/pythonTM

# link default output log to new log with subject identfier
ln -f ${dir_in}/output_logs/pythonTM_${SLURM_JOB_ID}.out ${dir_in}/output_logs/pythonTM_sub-${SUB}_ses-${SES}_task-${TASK}.out
ln -f ${dir_in}/output_logs/pythonTM_${SLURM_JOB_ID}.err ${dir_in}/output_logs/pythonTM_sub-${SUB}_ses-${SES}_task-${TASK}.err

# update output log with script details
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "job: ${SLURM_JOB_ID}"
echo "sub: sub-${SUB}"
echo "ses: ses-${SES}"
echo "task: task-${TASK}"
echo ""

# make subject output folder
if [ ! -d ${output_folder}/sub-${SUB}/ses-${SES} ]; then
    mkdir -p ${output_folder}/sub-${SUB}/ses-${SES}/func
fi

# user specified options
# dconn_infile: full path to input dconn based on command line arguements
# template_infile: full path to template which should match the units of input dconn (r values, fisher z, z-scored)
# template_minthreshold: applied to network template to select most strongly associated greyordinates
# refine_template_matching: if true, will refine template matching for motor_grays, dmn_grays and sal_grays
# refine_minthreshold: what refinement threshold will be applied to the motor, dmn, and sal grays
dconn_infile=${input_folder}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH.dconn.nii
template_infile=/home/faird/shared/code/internal/analytics/compare_matrices_python/code/template_matching_python/support_files/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat
template_minthreshold=1.00
# WARNING: REFINEMENT IS NOT RECOMMENDED
refine_template_matching=True
refine_minthreshold=3.0

# set output filepaths
# matfile_out is a matlab matrix which has:
# i.  sum of squares between greyordinates and all networks
# ii. network assignments following winner-take-all assessment of effect sizes
matfile_out=${output_folder}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.mat
# dscalar_out will be a cifti scalar to visualize winner-take-all network assignments
# this requires an input dscalar to use as a template (for header and for colors)
dscalar_template_infile=/home/faird/shared/code/internal/analytics/compare_matrices_python/code/template_matching_python/support_files/dscalar-template_data-mint_networks-15_colors-power_cifti-2.dscalar.nii
dscalar_outfile=${output_folder}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.dscalar.nii
dscalar_outfile_scan=${output_folder}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonScanTM.dscalar.nii

# path to the python code
code_path=/home/faird/shared/code/internal/analytics/compare_matrices_python/code/template_matching_python/template_matching_v6.0.py

# run template matching
source /home/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
python -u ${code_path} ${dconn_infile} ${template_infile} ${template_minthreshold} ${refine_template_matching} ${refine_minthreshold} ${matfile_out} ${dscalar_template_infile} ${dscalar_outfile} ${dscalar_outfile_scan}

echo "script finished at at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"

#remove original output logs
rm ${dir_in}/output_logs/pythonTM_${SLURM_JOB_ID}.out
rm ${dir_in}/output_logs/pythonTM_${SLURM_JOB_ID}.err