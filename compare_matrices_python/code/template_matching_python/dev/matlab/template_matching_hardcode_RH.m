%% Script overview 
% This is a hardcoded version of the original matlab template matching implimentation 
% script is potentially useful for comparision between matlab and python implimentations

%This code assigns each grayordinate a network based on the sptialally
%similiarty to a previously defined independent template. (Look at
%Make_Cifti_template_RH.m).  

%inputs are:
%dconn_filename is the path to your dconn
%data_type = "parcellated" or "dense" connectivity matrix
%template_path = path to .mat file that has th network templates.
%transform_data =  if you want to convert you data before comparing to your template, use can use 1 of 2 transformations: 'Convert_FisherZ_to_r' or 'Convert_r_to_Pearons' or 'Convert_to_Zscores' or use no tranformation.
%output_cifti_name = output_cifti_name pretty clear.
%cifti_output_folder = your project directory
%wb_command = path to run HCP workbench command.
%make_cifti_from_results = set to 1 if you want to save your results as a cifti. 0 will not save anything.
%allow_overlap = set to 1 if you're using overlapping networks in your cifti (Your input networks file will likely be a .dtseries.nii)
%overlap_method =  currently, the only supported method is "smooth_then_derivative"
%surface_only = set to 1 if you only want to generate assignments for the
%cortex (asssuming 59412 grayordinates.), otherwise set to zero and
%the output will be a dsclar of the standard size (91282).
%already_surface_only = set to 1 if your dconn is already 59412 x 59412.

% outputs are
% dscalar file of the assignments  and a cleaned file with the small islands removed.
% a dtseries file with the overlapping networks.
% a .mat file of the weights.

dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn';
cifti_output_folder='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/template_matching/derivatives';
dconn_filename='/scratch.global/kgodfrey/DCAN_WIP/dconnZ/sub-10938/ses-combined/func/sub-10938_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH.dconn.nii';
template_path='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/support_files/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat';
output_cifti_name='sub-10938_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm_zscoreRH_matlabTM_RH';
TEMPLATEMINIMUM = 1.00;
SCANTEMPLATEMINIMUM_thresholds = 3;
wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/bin_linux64/wb_command';
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'     ''    'SCAN'};

% parameters
% adding paths for this function
this_code = which('template_matching_RH');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files']; %find support files in the code directory.
addpath(genpath(support_folder));
settings=settings_comparematrices;%
np=size(settings.path,2);

disp('Attempting to add neccesaary paths and functions.')
addpath(genpath('/home/faird/shared/code/external/utilities/cifti-matlab'));
warning('on')
wb_command=settings.path_wb_c; %path to wb_command

disp('opening template...')
disp(['template minimum is set at ' num2str(TEMPLATEMINIMUM)]);
load(template_path,'seed_matrix'); % specify the seedmatrix variable to load just in case there is a variable with the same name.
cifti_template_mat_full =seed_matrix;
cifti_template_mat_full(cifti_template_mat_full<= TEMPLATEMINIMUM) = nan;

disp('opening subject dconn...')
subject_cii=ciftiopen(char(dconn_filename), wb_command); %dconn path
corr_mat_full = single(subject_cii.cdata);
clear subject_cii %save memory

disp('Calculating similarity (eta) to template')
%%% compute eta similarity value b/w each vertex and template %%%
if exist('eta_to_template_vox','var') ~=1
    eta_to_template_vox = single(zeros(size(corr_mat_full,1),length(network_names)));
    disp(['your template is size : ' num2str(size(cifti_template_mat_full,1)) ' by '  num2str(size(cifti_template_mat_full,2)) ]);
    disp(['your correlation matrix is size : ' num2str(size(corr_mat_full,1)) ' by '  num2str(size(corr_mat_full,2)) ]);
    for i=1:size(corr_mat_full,1)
        if rem(i,5000)==0
            disp([' Calculating voxel ' num2str(i)]);toc;
        end
        for j=1:length(network_names)
            if j==4 || j ==6 || j==17
                continue
            end
            %%% compute an eta value for each voxel for each network (from fran's etacorr script) %%%
            %goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(j,:)));
            goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(:,j))');
            cmap = corr_mat_full(i,goodvox)';
            %tmap = cifti_template_mat_full(j,goodvox)';
            tmap = cifti_template_mat_full(goodvox,j);
            Mgrand  = (mean(mean(tmap)) + mean(mean(cmap)))/2;
            Mwithin = (tmap+cmap)/2;
            SSwithin = sum(sum((tmap-Mwithin).*(tmap-Mwithin))) + sum(sum((cmap-Mwithin).*(cmap-Mwithin)));
            SStot    = sum(sum((tmap-Mgrand ).*(tmap-Mgrand ))) + sum(sum((cmap-Mgrand ).*(cmap-Mgrand )));
            eta_to_template_vox(i,j) = 1 - SSwithin/SStot;
            
            clear cmap tmap Mgrand Mwithin SSwithin SStot goodvox
        end
    end
end

%%% winner-take-all: highest eta value is network that voxel will be assigned to %%%
[~,new_subject_labels] = max(eta_to_template_vox,[],2); %find max for template matching
clear i

% increase the Z-score threshold to identify the scan network.
if size(eta_to_template_vox,2) == 18
    %make copies of the winner-take-all-results and put them
    %into an array:
    %new_subject_labels_scan_thresholds=repmat(new_subject_labels,size(SCANTEMPLATEMINIMUM_thresholds),1);
    new_subject_labels_scan_thresholds=new_subject_labels;
    
    disp(['testing scan threshold: ' num2str(SCANTEMPLATEMINIMUM_thresholds)]);
    motor_grays=find(new_subject_labels_scan_thresholds == 10 | new_subject_labels_scan_thresholds == 11 | new_subject_labels_scan_thresholds == 18); % find the networks that belong to SMd(10), SMl(11), or SCAN (18).
    
    % recalculate the seed map usng a higher threshold.
    cifti_template_mat_full_scan =seed_matrix;
    cifti_template_mat_full_scan(cifti_template_mat_full_scan<= SCANTEMPLATEMINIMUM_thresholds) = nan;
    
    disp(['Recalculating similarity (eta) to template using only motor networks: n= ' num2str(size(motor_grays,1)) ' grayordinates.'])
    %%% compute eta similarity value b/w each vertex and template %%%
    eta_to_template_vox_scan = single(zeros(size(motor_grays,1),length(network_names)));
    for i=1:size(eta_to_template_vox_scan,1)
        if rem(i,5000)==0
            disp([' Calculating voxel ' num2str(i)]);toc;
        end
        for j=1:length(network_names)
            if j==4 || j ==6 || j==17
                continue
            end
            %%% compute an eta value for each voxel for each network (from fran's etacorr script) %%%
            %goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(j,:)));
            goodvox = (~isnan(corr_mat_full(motor_grays(i),:)) & ~isnan(cifti_template_mat_full_scan(:,j))');
            cmap = corr_mat_full(motor_grays(i),goodvox)';
            %tmap = cifti_template_mat_full(j,goodvox)';
            tmap = cifti_template_mat_full_scan(goodvox,j);
            Mgrand  = (mean(mean(tmap)) + mean(mean(cmap)))/2;
            Mwithin = (tmap+cmap)/2;
            SSwithin = sum(sum((tmap-Mwithin).*(tmap-Mwithin))) + sum(sum((cmap-Mwithin).*(cmap-Mwithin)));
            SStot    = sum(sum((tmap-Mgrand ).*(tmap-Mgrand ))) + sum(sum((cmap-Mgrand ).*(cmap-Mgrand )));
            eta_to_template_vox_scan(i,j) = 1 - SSwithin/SStot;
            
            clear cmap tmap Mgrand Mwithin SSwithin SStot goodvox
        end
    end
    
    disp('Modifying motor grayordinates with the new higher threshold results...' )
    eta_to_template_vox_modified=eta_to_template_vox;
    eta_to_template_vox_modified(motor_grays,:) = eta_to_template_vox_scan;
    
    %%% winner-take-all: highest eta value is network that voxel will be assigned to %%%
    [~,new_subject_labels_scan_thresholds] = max(eta_to_template_vox_modified,[],2); %find max for template matching
end
clear goodvox i temp
            
disp(['Saving .mat file: ' cifti_output_folder '/' output_cifti_name '.mat'])
if size(eta_to_template_vox,2) == 18
    save([cifti_output_folder '/' output_cifti_name '.mat'],'eta_to_template_vox','new_subject_labels','network_names','eta_to_template_vox_modified','new_subject_labels_scan_thresholds','-v7.3')
else
    save([cifti_output_folder '/' output_cifti_name '.mat'],'eta_to_template_vox','new_subject_labels','network_names','-v7.3')
end

disp('saving no refinement dscalar to cifti')
dscalar_template_path='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/template_matching/support_files/sub-10227_ses-combined_task-cross_grims_recoloredCifti2.dscalar.nii';
saving_template =ciftiopen(dscalar_template_path, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
saving_template.cdata = single(new_subject_labels);
outfile_path = [cifti_output_folder '/' output_cifti_name '.dscalar.nii']
disp(outfile_path)
cifti_write(saving_template,[cifti_output_folder '/' output_cifti_name '.dscalar.nii'])
clear dscalar_template_path saving_template outfile_path

disp('saving refinement dscalar to cifti')
dscalar_template_path = '/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/template_matching/support_files/sub-10227_ses-combined_task-cross_grims_recoloredCifti2.dscalar.nii';
saving_template = ciftiopen(dscalar_template_path, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
saving_template.cdata = single(new_subject_labels_scan_thresholds);
outfile_path = [cifti_output_folder '/' output_cifti_name '_scanthresh3.dscalar.nii']
disp(outfile_path)
cifti_write(saving_template,outfile_path)
clear dscalar_template_path saving_template outfile_path