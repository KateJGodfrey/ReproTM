% code for simple template matching
% assumes that you have already created Z-scored dconns

dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn';
cifti_output_folder='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/v6_interpolated';
dconn_file='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_fisherZ_zscoreRH.dconn.nii';
template_path='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/support_files/seedmaps_subs_withsmoothed_dtseries_n141_all_networksZscored.mat';
output_cifti_name='sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_fisherZ_zscoreRH_matlabTM';
TEMPLATEMINIMUM = 1;
SCANTEMPLATEMINIMUM_thresholds = 3;
wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/bin_linux64/wb_command';
network_names = {   'DMN'    'Vis'    'FP'    ''    'DAN'     ''      'VAN'   'Sal'    'CO'    'SMd'    'SMl'    'Aud'    'Tpole'    'MTL'    'PMN'    'PON'     ''    'SCAN'};
allow_overlap = 0;

%% add paths to this function
this_code = which('template_matching_simple');
[code_dir,~] = fileparts(this_code);
support_folder=[code_dir '/support_files'];
addpath(genpath(support_folder));
% settings=settings_comparematrices;%
% np=size(settings.path,2);
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
addpath(genpath('/home/faird/shared/code/external/utilities/cifti-matlab'));

%% load the template
disp(['loading the following template: ',template_path])
disp(' ')
load(template_path);
% assign the seed matrix to a variable
cifti_template_mat_full = seed_matrix;
% exclude seeds below a threshold
disp(['template minimum is set at ' num2str(TEMPLATEMINIMUM)]);
disp('excluding seeds below template minimum')
disp(' ')
cifti_template_mat_full(cifti_template_mat_full<= TEMPLATEMINIMUM) = nan;


%% load the dconn
dconn_filename=strcat(dir_in,dconn_file);
subject_cii=ciftiopen(char(dconn_filename), wb_command); 
corr_mat_full = single(subject_cii.cdata);
% if range(corr_mat_full)>2
%     disp('The range of input cifti is greater than 2.  Your correlation matrix is probably Fisher Z tranformed (or Z-scored). Ensure that your template is tranformed similarly or set "Convert_FisherZ_to_r" to "1" to have it automatically tranformed to Pearson.');
% elseif range(corr_mat_full) <=2
%     disp('The range of input cifti is less than (or equal to) 2.  Your correlation matrix is probably in Pearson Correlation. Ensure that your template is also a pearson correlation.')
% end

clear subject_cii % save memory

disp('Calculating similarity (eta) to template')
%%% compute eta similarity value b/w each vertex and template %%%
eta_to_template_vox = single(zeros(size(corr_mat_full,1),length(network_names)));
for i=1:size(corr_mat_full,1)
    if rem(i,5000)==0
        disp([' Calculating voxel ' num2str(i)]);toc;
    end
    for j=1:length(network_names)
        if j==4 || j ==6
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
        eta = 1 - SSwithin/SStot;
        eta_to_template_vox(i,j) = 1 - SSwithin/SStot;
                
        clear cmap tmap Mgrand Mwithin SSwithin SStot goodvox
    end
end

%%% winner-take-all: highest eta value is network that voxel will be assigned to %%%
%new_subject_labels=strcat(dir_in,'/',output_cifti_name)
[~, new_subject_labels] = max(eta_to_template_vox,[],2); %find max for template matching
clear i

%% increase Z-score threshold to identify the scan network
% increase the Z-score threshold to identify the scan network.
if size(eta_to_template_vox,2) == 18
    %make copies of the winner-take-all-results and put them
    %into an array:
    new_subject_labels_scan_thresholds=repmat(new_subject_labels,size(SCANTEMPLATEMINIMUM_thresholds),1);
    
    for t=1:size(SCANTEMPLATEMINIMUM_thresholds,2)
        disp(['testing scan threshold: ' num2str(SCANTEMPLATEMINIMUM_thresholds(t))]);
        motor_grays=find(new_subject_labels_scan_thresholds(:,t) == 10 | new_subject_labels_scan_thresholds(:,t) == 11 | new_subject_labels_scan_thresholds(:,t) == 18); % find the networks that belong to SMd(10), SMl(11), or SCAN (18).
        
        % recalculate the seed map usng a higher threshold.
        cifti_template_mat_full_scan =seed_matrix;
        cifti_template_mat_full_scan(cifti_template_mat_full_scan<= SCANTEMPLATEMINIMUM_thresholds(t)) = nan;
        
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
    end
    eta_to_template_vox_modified=eta_to_template_vox;
    eta_to_template_vox_modified(motor_grays,:) = eta_to_template_vox_scan;

    %DMN refinement
    for t=1:size(SCANTEMPLATEMINIMUM_thresholds,2)
        disp(['testing scan threshold: ' num2str(SCANTEMPLATEMINIMUM_thresholds(t))]);
        DMN_grays=find(new_subject_labels_scan_thresholds(:,t) == 1 | new_subject_labels_scan_thresholds(:,t) == 15 | new_subject_labels_scan_thresholds(:,t) == 16); % find the networks that belong to SMd(10), SMl(11), or SCAN (18).
        
        % recalculate the seed map usng a higher threshold.
        cifti_template_mat_full_scan =seed_matrix;
        cifti_template_mat_full_scan(cifti_template_mat_full_scan<= SCANTEMPLATEMINIMUM_thresholds(t)) = nan;
        
        disp(['Recalculating similarity (eta) to template using only motor networks: n= ' num2str(size(DMN_grays,1)) ' grayordinates.'])
        %%% compute eta similarity value b/w each vertex and template %%%
        eta_to_template_vox_DMN = single(zeros(size(DMN_grays,1),length(network_names)));
        for i=1:size(eta_to_template_vox_DMN,1)
            if rem(i,5000)==0
                disp([' Calculating voxel ' num2str(i)]);toc;
            end
            for j=1:length(network_names)
                if j==4 || j ==6 || j==17
                    continue
                end
                %%% compute an eta value for each voxel for each network (from fran's etacorr script) %%%
                %goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(j,:)));
                goodvox = (~isnan(corr_mat_full(DMN_grays(i),:)) & ~isnan(cifti_template_mat_full_scan(:,j))');
                cmap = corr_mat_full(DMN_grays(i),goodvox)';
                %tmap = cifti_template_mat_full(j,goodvox)';
                tmap = cifti_template_mat_full_scan(goodvox,j);
                Mgrand  = (mean(mean(tmap)) + mean(mean(cmap)))/2;
                Mwithin = (tmap+cmap)/2;
                SSwithin = sum(sum((tmap-Mwithin).*(tmap-Mwithin))) + sum(sum((cmap-Mwithin).*(cmap-Mwithin)));
                SStot    = sum(sum((tmap-Mgrand ).*(tmap-Mgrand ))) + sum(sum((cmap-Mgrand ).*(cmap-Mgrand )));
                eta_to_template_vox_DMN(i,j) = 1 - SSwithin/SStot;
                
                clear cmap tmap Mgrand Mwithin SSwithin SStot goodvox
            end
        end
        
        disp('Modifying motor grayordinates with the new higher threshold results...' )
        %eta_to_template_vox_modified=eta_to_template_vox;
        eta_to_template_vox_modified(DMN_grays,:) = eta_to_template_vox_DMN;
        
        
        %refine salience
        for t=1:size(SCANTEMPLATEMINIMUM_thresholds,2)
            disp(['testing scan threshold: ' num2str(SCANTEMPLATEMINIMUM_thresholds(t))]);
            sal_grays=find(new_subject_labels_scan_thresholds(:,t) == 8 | new_subject_labels_scan_thresholds(:,t) == 9); % find the networks that belong to SMd(10), SMl(11), or SCAN (18).
            
            % recalculate the seed map usng a higher threshold.
            cifti_template_mat_full_scan =seed_matrix;
            cifti_template_mat_full_scan(cifti_template_mat_full_scan<= SCANTEMPLATEMINIMUM_thresholds(t)) = nan;
            
            disp(['Recalculating similarity (eta) to template using only motor networks: n= ' num2str(size(sal_grays,1)) ' grayordinates.'])
            %%% compute eta similarity value b/w each vertex and template %%%
            eta_to_template_vox_sal = single(zeros(size(sal_grays,1),length(network_names)));
            for i=1:size(eta_to_template_vox_sal,1)
                if rem(i,5000)==0
                    disp([' Calculating voxel ' num2str(i)]);toc;
                end
                for j=1:length(network_names)
                    if j==4 || j ==6 || j==17
                        continue
                    end
                    %%% compute an eta value for each voxel for each network (from fran's etacorr script) %%%
                    %goodvox = (~isnan(corr_mat_full(i,:)) & ~isnan(cifti_template_mat_full(j,:)));
                    goodvox = (~isnan(corr_mat_full(sal_grays(i),:)) & ~isnan(cifti_template_mat_full_scan(:,j))');
                    cmap = corr_mat_full(sal_grays(i),goodvox)';
                    %tmap = cifti_template_mat_full(j,goodvox)';
                    tmap = cifti_template_mat_full_scan(goodvox,j);
                    Mgrand  = (mean(mean(tmap)) + mean(mean(cmap)))/2;
                    Mwithin = (tmap+cmap)/2;
                    SSwithin = sum(sum((tmap-Mwithin).*(tmap-Mwithin))) + sum(sum((cmap-Mwithin).*(cmap-Mwithin)));
                    SStot    = sum(sum((tmap-Mgrand ).*(tmap-Mgrand ))) + sum(sum((cmap-Mgrand ).*(cmap-Mgrand )));
                    eta_to_template_vox_sal(i,j) = 1 - SSwithin/SStot;
                    
                    clear cmap tmap Mgrand Mwithin SSwithin SStot goodvox
                end
            end
            
            disp('Modifying motor grayordinates with the new higher threshold results...' )
            %eta_to_template_vox_modified=eta_to_template_vox;
            eta_to_template_vox_modified(sal_grays,:) = eta_to_template_vox_sal;
        end
        %%% winner-take-all: highest eta value is network that voxel will be assigned to %%%
        %[~, new_subject_labels] = max(eta_to_template_vox_modified,[],2); %find max for template matching
        [~, new_subject_labels_scan_thresholds(:,t)] = max(eta_to_template_vox_modified,[],2);
    end

    clear goodvox i temp

    if exist('allow_overlap','var') == 1 && allow_overlap == 1
        if allow_overlap == 1
            disp('Calculating overlap')
            MuI_threshhold_all_networks = findoverlapthreshold(eta_to_template_vox,network_names,Zscore_eta, overlap_method);
        else
        end
    else
    end

    toc
    disp(['Saving .mat file: ' cifti_output_folder '/' output_cifti_name '.mat'])
    if size(eta_to_template_vox,2) == 18
        save([cifti_output_folder '/' output_cifti_name '.mat'],'eta_to_template_vox','eta_to_template_vox_modified','new_subject_labels','network_names','new_subject_labels_scan_thresholds','-v7.3')
    else
        save([cifti_output_folder '/' output_cifti_name '.mat'],'eta_to_template_vox','new_subject_labels','network_names','-v7.3') 
    end

    % switch transform_data
    %     case 'Convert_to_Zscores'
    %         unix(['rm -f ' char(Zdconn) ])
    %     otherwise
    % end

    % disp('saving file to cifti')
    % if surface_only ==1
    %     saving_template = ciftiopen(settings.path{12}, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
    % else %assume 91282
    %     saving_template = ciftiopen(settings.path{8}, wb_command); % don't forget to load in a gifti object, or  else saving_template will be interpreted as a struct.
    % end
    
    % if  size(eta_to_template_vox,2) == 18
    %     new_subject_labels = new_subject_labels_scan_thresholds;
    % end

end
