%function[output_cifti_name] = Zscore_dconn(input_cifti_name,output_cifti_name,path_wb_c)

dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn'
input_cifti_name=strcat(dir_in,'/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_fisherZ.dconn.nii')
output_cifti_name=strcat(dir_in, '/sub-10227/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_fisherZ.dconn.nii')

path_wb_c='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/bin_linux64/wb_command';
%%addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/Matlab_CIFTI'))

%NOTE: Cannot use the following path because it loads cifti data as a struct rather than a gifti object.
%addpath(genpath('/home/faird/shared/code/external/utilities/cifti-matlab/')) % CHANGE MADE BY CristianM 03/15/2021
amIdeployed = isdeployed();
amIdeployed=num2str(double(amIdeployed));
disp(['is deployed equals: ' amIdeployed]);


addpath(genpath('/home/faird/shared/code/internal/utilities/Matlab_CIFTI/')) %-CHANGE MADE BY RobertH 04/30/2021
%addpath(genpath('/mnt/max/shared/code/internal/utilities/CIFTI/'))
%%addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/external/utilities/gifti-1.6'))
addpath(genpath('/home/faird/shared/code/external/utilities/gifti-1.6/gifti-1.6/'))  % CHANGE MADE BY CristianM 03/15/2021

%path_wb_c='LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command';
%path_wb_c='/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-9253ac2/bin_rh_linux64/wb_command';

%path_wb_c='/home/exacloud/lustre1/fnl_lab/code/external/utilities/workbench-1.2.3-HCP/bin_rh_linux64/wb_command';
%path_wb_c='/panfs/roc/msisoft/workbench/1.4.2/bin/wb_command';% CHANGE MADE BY CristianM 03/15/2021
%path_wb_c='/home/faird/shared/code/external/utilities/workbench/1.4.2/workbench/bin_rh_linux64/wb_command'; %-CHANGE MADE BY RobertH 04/30/2021
if strcmp(output_cifti_name,'inferred') == 1
    short_file_name = char(input_cifti_name(1:end-10)); %(i.e. remove ".dconn.nii");
    output_cifti_name = [short_file_name 'Zscored'];
else
    disp('output_cifti_name must be provided or type "inferred" and it will be placed in the same location.');
end

if exist([output_cifti_name '.dconn.nii'],'file') ~=0
    disp(['Zscored file already found at this location: '  output_cifti_name '.dconn.nii']);
    disp('No  need to remake.');
    output_cifti_name = [output_cifti_name '.dconn.nii']; %return the output name with the extension.    
else
    %addpath(genpath('/mnt/max/shared/code/external/utilities/HCP_CCA'))
    %wb_command='wb_command';
    %subjectlist=importdata('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/palm_ADHD_ptseries.conc_pconn_of_ptseries_5_minutes_of_data_at_FD_0.2.conc');
    %Controlsubjectlist=importdata('/mnt/max/shared/code/internal/utilities/hcp_comm_det_damien/palm_Control_ptseries.conc_pconn_of_ptseries_5_minutes_of_data_at_FD_0.2.conc');
    disp('loading cifti');
    tic
    newcii = ciftiopen(input_cifti_name,path_wb_c);
    dconn=single(newcii.cdata);
    disp('first few elements of original dconn');
    dconn(1:5,1:5)
    toc
    %clear newcii; %save memory
    %     for i = 1:length(subjectlist) %check for subjects
    %         if exist(subjectlist{i}) == 0
    %             disp(['Subject series ' num2str(i) ' does not exist'])
    %             %return
    %         else
    %         end
    %     end
    %A=size(subjectlist,1);
    if exist('exclude_lower_triangle','var') == 1 % not debugged
        LL_mat_size = 29696; %hardcode - number of grayordinates
        LR_mat_size = 29716;
        SS_mat_size = 31870;
        
        LLlowern=mat_size-1; %used for lower calculation of triangle
        LLGausse_num=((lowern*(lowern+1))/2);
        LLsubjectmat=zeros(LLGausse_num,1);
        
        LRlowern=mat_size-1; %used for lower calculation of triangle
        LRGausse_num=((lowern*(lowern+1))/2);
        LRsubjectmat=zeros(LLGausse_num,1);
        
        SSlowern=mat_size-1; %used for lower calculation of triangle
        SSGausse_num=((lowern*(lowern+1))/2);
        SSsubjectmat=zeros(LLGausse_num,1);
        
        for i=1:size(subjectlist) % go through each subject
            subjectCIFTI= ciftiopen(char(subjectlist(i)),wb_command);
            %subj_conn_column=zeros((sz*(sz+1)/2),1); %use the Gausse's formula (n*(n+1)/2) to make the appropraite number of columns
            %for j=1:size(sz)
            %    subj_conn_column((1:sz)*j,1)=(subjectCIFTI.cdata(j+1:sz,j));
            %end
            n = size(subjectCIFTI.cdata,1);
            v = subjectCIFTI.cdata(find(tril(ones(n,n),-1)));
            ADHD_subjectmat(:,i)=v;
        end
        
    end
    
    disp ('sectioning dconn')
    LL = dconn(1:29696,1:29696);
    LR = dconn(29697:59412,1:29696);
    LS = dconn(59413:91282,1:29696);
    RL = dconn(1:29696,29697:59412);
    RR = dconn(29697:59412,29697:59412);
    RS = dconn(59413:91282,29697:59412);
    SL = dconn(1:29696,59413:91282);
    SR = dconn(29697:59412,59413:91282);
    SS = dconn(59413:91282,59413:91282);
    
    disp('resphaping nonants to calcuate Z scores')
    ZLL = zscore(reshape(LL,size(LL,1)*size(LL,2),1));
    ZLLmat = reshape(ZLL,size(LL,1),size(LL,2)); clear ZLL LL
    disp('ZLL');
    
    ZLR = zscore(reshape(LR,size(LR,1)*size(LR,2),1));
    ZLRmat = reshape(ZLR,size(LR,1),size(LR,2)); clear ZLR LR
    disp('ZLR');
    
    ZLS = zscore(reshape(LS,size(LS,1)*size(LS,2),1));
    ZLSmat = reshape(ZLS,size(LS,1),size(LS,2)); clear ZLS LS
    disp('ZLS');
    
    ZRL = zscore(reshape(RL,size(RL,1)*size(RL,2),1));
    ZRLmat = reshape(ZRL,size(RL,1),size(RL,2)); clear ZRL RL
    disp('ZRL');
    
    ZRR = zscore(reshape(RR,size(RR,1)*size(RR,2),1));
    ZRRmat = reshape(ZRR,size(RR,1),size(RR,2)); clear ZRR RR
    disp('ZRR');
    
    ZRS = zscore(reshape(RS,size(RS,1)*size(RS,2),1));
    ZRSmat = reshape(ZRS,size(RS,1),size(RS,2)); clear ZRS RS
    disp('ZRS');
    
    ZSL = zscore(reshape(SL,size(SL,1)*size(SL,2),1));
    ZSLmat = reshape(ZSL,size(SL,1),size(SL,2)); clear ZSL SL
    disp('ZSL');
    
    ZSR = zscore(reshape(SR,size(SR,1)*size(SR,2),1));
    ZSRmat = reshape(ZSR,size(SR,1),size(SR,2)); clear ZSR SR
    disp('ZSR');
    
    ZSS = zscore(reshape(SS,size(SS,1)*size(SS,2),1));
    ZSSmat = reshape(ZSS,size(SS,1),size(SS,2)); clear ZSS SS
    disp('ZSS');
    
    disp('rewriting matrix')
    newdconn=single(zeros(91282,91282));
    newdconn(1:29696,1:29696) = ZLLmat;
    newdconn(29697:59412,1:29696) = ZLRmat;
    newdconn(59413:91282,1:29696) = ZLSmat;
    newdconn(1:29696,29697:59412) = ZRLmat;
    newdconn(29697:59412,29697:59412) = ZRRmat;
    newdconn(59413:91282,29697:59412) = ZRSmat;
    newdconn(1:29696,59413:91282) = ZSLmat;
    newdconn(29697:59412,59413:91282) = ZSRmat;
    newdconn(59413:91282,59413:91282) = ZSSmat;
    
    disp('first few rows and columns of new dconn:')
    newdconn(1:5,1:5)
    
    clear dconn
    %newcii = ciftiopen('/mnt/max/shared/projects/midnight_scan_club/info_map/Results/Community_Detection_Min_Dist_30_TieDen_0.03_MinNet_Size_400_MinReg_Size_30/MSC_Cesna_test_dthalf1.conc_all_frames_at_FD_thresh_0.2_and_smoothing_2.55_AVG.dconn.nii/dconn/MSC_Cesna_test_dthalf1.conc_all_frames_at_FD_thresh_0.2_and_smoothing_2.55_AVG.dconn.nii','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/bin/wb_command');
    newcii.cdata = newdconn;
    %addpath('/mnt/max/shared/code/internal/utilities/corr_pt_dt/support_files');
    disp('Saving new Zscored dconn')
    % save(newcii, [output_cifti_name '.gii'], 'ExternalFileBinary')
    % disp('Converting Zscored dconn .gii to .nii')
    % unix([path_wb_c ' -cifti-convert -from-gifti-ext ' output_cifti_name '.gii ' output_cifti_name '.dconn.nii' ]);
    % disp('Removing .gii')
    % unix(['rm -f ' output_cifti_name '.gii']);
    % unix(['rm -f ' output_cifti_name '.dat']);
    
    
    
    
    save(newcii, [output_cifti_name '.gii'], 'ExternalFileBinary')
    disp('Converting Zscored dconn .gii to .nii')
    unix([path_wb_c ' -cifti-convert -from-gifti-ext ' output_cifti_name '.gii ' output_cifti_name '.dconn.nii' ]);
    disp('Removing .gii')
    unix(['rm -f ' output_cifti_name '.gii']);
    unix(['rm -f ' output_cifti_name '.dat']);
    
    output_cifti_name = [output_cifti_name '.dconn.nii']; %return the output name with the extension.
end

% if exist('plot_distributions','var') == 1 % not debugged
    
%     figure()
%     subplot(3,3,1) %FD plot
%     histogram(LL,'BinLimits',[-0.3,0.3],'facecolor',[1 0 0], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('LL');
%     subplot(3,3,2) %FD plot
%     histogram(LR,'BinLimits',[-0.3,0.3],'facecolor',[0 1 0], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('LR');
%     subplot(3,3,3) %FD plot
%     histogram(LS,'BinLimits',[-0.3,0.3],'facecolor',[1 0.5 0], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('LS');
%     subplot(3,3,4) %FD plot
%     histogram(RL,'BinLimits',[-0.3,0.3],'facecolor',[0 1 0], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('RL');
%     subplot(3,3,5) %FD plot
%     histogram(RR,'BinLimits',[-0.3,0.3],'facecolor',[0 1 1], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('RR');
%     subplot(3,3,6) %FD plot
%     histogram(RS,'BinLimits',[-0.3,0.3],'facecolor',[0.5 0 0.5], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('RS');
%     subplot(3,3,7) %FD plot
%     histogram(SL,'BinLimits',[-0.3,0.3],'facecolor',[1 0.5 0], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('SL');
%     subplot(3,3,8) %FD plot
%     histogram(SR,'BinLimits',[-0.3,0.3],'facecolor',[0.5 0 0.5], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('SR');
%     subplot(3,3,9) %FD plot
%     histogram(SS,'BinLimits',[-0.3,0.3],'facecolor',[0 0 1], 'facealpha', 0.5, 'edgecolor', 'none')
%     title('SS');
    
%     figure()
%     %         histogram(LL,'BinLimits',[-0.3,0.3],'facecolor',[1 0 0], 'facealpha', 0.5, 'edgecolor', 'none'); hold on
%     %         histogram(LR,'BinLimits',[-0.3,0.3],'facecolor',[0 1 0], 'facealpha', 0.5, 'edgecolor', 'none'); hold on
%     %         histogram(SS,'BinLimits',[-0.3,0.3],'facecolor',[0 0 1], 'facealpha', 0.5, 'edgecolor', 'none'); hold on
    
%     histogram([nonzeros(triu(LL,1)); nonzeros(triu(RR,1)); reshape(LR,size(LR,1)*size(LR,2),1) ],'BinLimits',[-0.3,0.3],'facecolor',[1 0 0], 'facealpha', 0.5, 'edgecolor', 'none'); hold on
%     histogram([reshape(RS,size(RS,1)*size(RS,2),1); reshape(LS,size(LS,1)*size(LS,2),1)],'BinLimits',[-0.3,0.3],'facecolor',[0 1 0], 'facealpha', 0.5, 'edgecolor', 'none'); hold on
%     histogram(SS,'BinLimits',[-0.3,0.3],'facecolor',[0 0 1], 'facealpha', 0.5, 'edgecolor', 'none'); hold on
    
%     legend({'cortical','cort2sub','subcortical'});
    
%     mat_size = 91282; %hardcode - number of parcellations
%     %dconns = single(zeros(mat_size,mat_size)); %build a 352 x 352 x subject matrix
%     %A=size(subjectlist,1);
%     A =1;
%     subjectlist = 1;
%     %used for unwrapping matrix
%     lowern=mat_size-1; %used for lower calculation of triangle
%     Gausse_num=((lowern*(lowern+1))/2);
%     %subjectmat=single(zeros(Gausse_num,A));
    
%     %for i=1:size(subjectlist) % go through each subject
%     %subjectCIFTI= ciftiopen(char(subjectlist(i)),wb_command);
%     %subj_conn_column=zeros((sz*(sz+1)/2),1); %use the Gausse's formula (n*(n+1)/2) to make the appropraite number of columns
%     %for j=1:size(sz)
%     %    subj_conn_column((1:sz)*j,1)=(subjectCIFTI.cdata(j+1:sz,j));
%     %end
%     %n = size(dconn,1);
%     %v = dconn(find(tril(ones(n,n),-1)));
%     %subjectmat(:,i)=v;
%     %dconns(:,:,i)=subjectCIFTI.cdata;
%     %disp(i);
%     %end
    
% end

% disp('Done coverting Dconn to Zscored Dconn.')

% end

