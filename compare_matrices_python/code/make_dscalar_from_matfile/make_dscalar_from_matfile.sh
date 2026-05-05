#!/bin/bash -l

# user specified inputs
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo " "

dir_in=/home/faird/shared/code/internal/analytics/compare_matrices_python
mat_filepath=${dir_in}/derivatives/pythonTM/sub-10938/ses-combined/func/sub-10938_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.mat
dscalar_outfile=${dir_in}/derivatives/pythonTM/sub-10938/ses-combined/func/sub-10938_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM_rvalues.dscalar.nii
dscalar_template_infile=/home/faird/shared/code/internal/analytics/compare_matrices_python/code/template_matching_python/support_files/dscalar-template_data-mint_networks-15_colors-power_cifti-2.dscalar.nii

# make dscalar using pearson correlation matching
code_path=/home/faird/shared/code/internal/analytics/compare_matrices_python/code/make_dscalar_from_matfile/make_dscalar_rvals.py

# make dscalar using sum of squares matching
code_path=/home/faird/shared/code/internal/analytics/compare_matrices_python/code/make_dscalar_from_matfile/make_dscalar.py

python -u ${code_path} ${mat_filepath} ${dscalar_outfile} ${dscalar_template_infile}

echo " " 
echo "script completed at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"