#!/bin/bash -l

#SBATCH -J clean_dscalar
#SBATCH --ntasks=4
#SBATCH --tmp=5gb
#SBATCH --mem=5gb
#SBATCH -t 0:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pandh015@umn.edu
#SBATCH -p msismall
#SBATCH -o /projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/output_logs/clean_dscalar_%j.out
#SBATCH -e /projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/output_logs/clean_dscalar_%j.err
#SBATCH -A midb_abcd

# user specified inputs from terminal
SUB=0A4P0LWM
SES=00A
TASK=rest

# this script runs a copy of the original matlab code
# which will clean up networks below a user-specified minimum size
# 	source code: /home/faird/shared/code/internal/analytics/compare_matrices_to_assign_networks/clean_dscalars_by_size.m
# this script takes as inputs:
# 	i.  template matched dscalar
# 	ii. user-specified minimum network size
# this script will output:
# 	i.  dscalar in the same folder with *_recolored_minsize

# update output log with script details
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "job: ${SLURM_JOB_ID}"
echo "sub: sub-${SUB}"
echo "ses: ses-${SES}"
echo "task: task-${TASK}"

# what directory are we in?
dir_in=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python
# what is the path to your input dscalar file for cleaning
dscalar_infile=${dir_in}/derivatives/pythonTM/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii_10_minutes_of_data_at_FD_0.2_zscoreRH_pythonTM.dscalar.nii
# what do you want to set as the minimum size for a region?
minsize=30
# what version of workbench do you want to use?
module load workbench/1.4.2

matlab_command="addpath(genpath('/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/clean_dscalars'));clean_dscalars_by_size_simple_v2('${dscalar_infile}',[],[],[],[],${minsize},[],0,0,1);"
matlab_exec="/common/software/install/migrated/matlab/R2019a/bin/matlab"

${matlab_exec} -nodisplay -nosplash -batch ${matlab_command}


# echo "script finished at: $(date +"%Y-%m-%d") $(date +"%r")"
# echo "****************************************************"
