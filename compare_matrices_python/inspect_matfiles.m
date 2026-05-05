%% overview
% helper script for use with matlab gui or command line to inspect your outputs
% useful to compare template matching outputs between matlab and python implementations

%% load data

% set working directory and add path
dir_in='/panfs/jay/groups/6/faird/shared/code/internal/analytics/compare_matrices_python';
addpath(genpath(dir_in));
% load matfiles from matlab and python template matching
matlabTM_matfile=load(strcat(dir_in,'/derivatives/pythonTM/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_matlabTM.mat'));
pythonTM_matfile=load(strcat(dir_in,'/derivatives/pythonTM/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.mat'));

%% compare eta matrices

% extract eta matrices from matfiles
matlab_etamatrix=matlabTM_matfile.eta_to_template_vox;
python_etamatrix=pythonTM_matfile.eta_to_template_vox;
python_rmatrix=pythonTM_matfile.r_to_template_vox;
% test similarity between matlab and python derived eta matrices (correlation)
corr(matlab_etamatrix(:),python_etamatrix(:),'rows','complete')

% extract refined eta matrices from matfiles 
matlab_etamatrix_refined=matlabTM_matfile.eta_to_template_vox_modified;
python_etamatrix_refined=pythonTM_matfile.eta_to_template_vox_modified;
% test similarity between matlab and python derived refined eta matrices (correlation)
corr(matlab_etamatrix_refined(:),python_etamatrix_refined(:),'rows','complete');

%% compare network assignments

% extract network assignment dscalars from matfiles
matlab_networks=matlabTM_matfile.new_subject_labels;
python_networks=pythonTM_matfile.new_subject_labels;
python_rnetworks=pythonTM_matfile.new_subject_labels_r;
% test similarity between matlab and python dscalars (# greys with differing assignments)
sum(matlab_networks ~= transpose(python_networks))
% test similarity between SS and Pearson r dscalars (# greys with differing assignments)
sum(transpose(python_rnetworks) ~= transpose(python_networks))
