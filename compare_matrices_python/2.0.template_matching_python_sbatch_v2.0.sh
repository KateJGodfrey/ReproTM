#!/bin/bash -l

#SBATCH -J pythonTM
#SBATCH --ntasks=2
#SBATCH --tmp=120gb
#SBATCH --mem=120gb
#SBATCH -t 2:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pandh015@umn.edu
#SBATCH -p msibigmem
#SBATCH -o /projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/output_logs/pythonTM_%j.out
#SBATCH -e /projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/output_logs/pythonTM_%j.err
#SBATCH -A midb_abcd

# update log
# v2.0.
# major updates:
# i. choice to threshold the template is now a user-specified option for use with template_matching_v6.0.py
# minor updates:
# i. updated user specified variable dscalar_outfile_scan to dscalar_outfile_refined

## TO DO LIST:
# i.   Try out with different templates for each sessions (00A, 02A, etc) i.e different age groups
# ii.  Test out other assignment methods ( consine similarity rather than sum of squares etc.)

SUB=0A4P0LWM
SES=00A
TASK=rest

# example command:
# sbatch ./2.0.template_matching_python_sbatch_v2.0.sh 10227 combined cross

# what directory are you in?
dir_in=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python

# what directory are your input dconns in?
input_folder='/scratch.global/pandh015/precision_maps_via_template_matching/dconnZ'
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
# refine_template_matching: if true, will refine template matching for somatomotor lateral, somatomotor dorsal, and scan
# refine_minthreshold: what template threshold will be applied during refinement of somatomotor and scan greys
dconn_infile=${input_folder}/sub-${SUB}/ses-${SES}/sub-${SUB}_ses-${SES}_task-${TASK}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii_10_minutes_of_data_at_FD_0.2_Zscore.dconn.nii
template_infile=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/template_matching_python/support_files/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat 
template_thresholding=yes      # options: yes, no 
template_minthreshold=1.00
refine_template_matching=yes    
refine_minthreshold=3.00        

# set output filepaths
# matfile_out is a matlab matrix which has:
# i.  sum of squares between greyordinates and all networks
# ii. network assignments following winner-take-all assessment of effect sizes
matfile_out=${output_folder}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii_10_minutes_of_data_at_FD_0.2_zscoreRH_pythonTM.mat
# dscalar_out will be a cifti scalar to visualize winner-take-all network assignments
# this requires an input dscalar to use as a template (for header and for colors)
dscalar_template_infile=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/template_matching_python/support_files/dscalar-template_data-mint_networks-15_colors-power_cifti-2.dscalar.nii
dscalar_outfile=${output_folder}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii_10_minutes_of_data_at_FD_0.2_zscoreRH_pythonTM.dscalar.nii
dscalar_outfile_refined=${output_folder}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii_10_minutes_of_data_at_FD_0.2_zscoreRH_pythonScanTM.dscalar.nii

# path to the python code
code_path=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/template_matching_python/template_matching_v6.0.py

# run template matching
source /home/faird/shared/code/external/envs/miniconda3/load_miniconda3.sh
python -u ${code_path} ${dconn_infile} ${template_infile} ${template_thresholding} ${template_minthreshold} ${refine_template_matching} ${refine_minthreshold} ${matfile_out} ${dscalar_template_infile} ${dscalar_outfile} ${dscalar_outfile_refined}

echo "script finished at at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"

#remove original output logs
rm ${dir_in}/output_logs/pythonTM_${SLURM_JOB_ID}.out
rm ${dir_in}/output_logs/pythonTM_${SLURM_JOB_ID}.err