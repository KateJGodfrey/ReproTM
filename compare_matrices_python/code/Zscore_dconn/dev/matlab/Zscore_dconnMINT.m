
module load matlab/R2019a
matlab

dir_in='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/dconn';
cifti_output_folder='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching';
dconn_file='/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_smoothed2mm_censoredFD2mm.dconn.nii';
output_cifti_name='sub-10227_ses-combined_task-taskMENORDICtrimmed_TM';
wb_command='/home/faird/shared/code/external/utilities/workbench/1.5.0/workbench/bin_linux64/wb_command';

%% add paths
addpath(genpath('/home/exacloud/lustre1/fnl_lab/code/internal/utilities/plotting-tools'));
addpath(genpath('/home/faird/shared/code/external/utilities/cifti-matlab'));

%% load the dconn
dconn_filename=strcat(dir_in,dconn_file);
subject_cii=ciftiopen(char(dconn_filename), wb_command); 
dconn = double(subject_cii.cdata);
%clear subject_cii  % save memory

%% can do some calculations on the dconn
% format long                      % format long prints in same notation as python
% matrix_sum=sum(nansum(dconn));   % sum of all the matrix elements
% num_nans=sum(sum(isnan(dconn))); % number of NaN values in matrix
    
disp ('sectioning dconn');
LL = dconn(1:29696,1:29696);
% LR = dconn(29697:59412,1:29696);
% LS = dconn(59413:91282,1:29696);
% RL = dconn(1:29696,29697:59412);
% RR = dconn(29697:59412,29697:59412);
% RS = dconn(59413:91282,29697:59412);
% SL = dconn(1:29696,59413:91282);
% SR = dconn(29697:59412,59413:91282);
% SS = dconn(59413:91282,59413:91282);
disp ('sectioning dconn done');

% can do some calculations on the sub matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_matrix = SS;
% matrix_sum=sum(nansum(test_matrix));   % sum of all the matrix elements
% num_nans=sum(sum(isnan(test_matrix))); % number of NaN values in matrix
% disp(strcat("number of NaN: ",num2str(num_nans)));
% disp(strcat("matrix sum: ",num2str(matrix_sum)));
% disp(strcat("matrix size: ",num2str(numel(test_matrix))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('resphaping nonants to calcuate Z scores')
ZLL = zscore(reshape(LL,size(LL,1)*size(LL,2),1));
ZLLmat = reshape(ZLL,size(LL,1),size(LL,2)); clear ZLL LL
disp('ZLL');

ZLR = zscore(reshape(LR,size(LR,1)*size(LR,2),1));
ZLRmat = reshape(ZLR,size(LR,1),size(LR,2)); %clear ZLR LR
disp('ZLR');

ZLS = zscore(reshape(LS,size(LS,1)*size(LS,2),1));
ZLSmat = reshape(ZLS,size(LS,1),size(LS,2)); %clear ZLS LS
disp('ZLS');

ZRL = zscore(reshape(RL,size(RL,1)*size(RL,2),1));
ZRLmat = reshape(ZRL,size(RL,1),size(RL,2)); %clear ZRL RL
disp('ZRL');

ZRR = zscore(reshape(RR,size(RR,1)*size(RR,2),1));
ZRRmat = reshape(ZRR,size(RR,1),size(RR,2)); %clear ZRR RR
disp('ZRR');

ZRS = zscore(reshape(RS,size(RS,1)*size(RS,2),1));
ZRSmat = reshape(ZRS,size(RS,1),size(RS,2)); %clear ZRS RS
disp('ZRS');

ZSL = zscore(reshape(SL,size(SL,1)*size(SL,2),1));
ZSLmat = reshape(ZSL,size(SL,1),size(SL,2)); %clear ZSL SL
disp('ZSL');

ZSR = zscore(reshape(SR,size(SR,1)*size(SR,2),1));
ZSRmat = reshape(ZSR,size(SR,1),size(SR,2)); %clear ZSR SR
disp('ZSR');

ZSS = zscore(reshape(SS,size(SS,1)*size(SS,2),1));
ZSSmat = reshape(ZSS,size(SS,1),size(SS,2)); %clear ZSS SS
disp('ZSS');

%% Can display information about matrices generated above %%%%%%%%%%%%%%
test_matrix = ZLLmat;
matrix_sum=sum(nansum(test_matrix));   % sum of all the matrix elements
num_nans=sum(sum(isnan(test_matrix))); % number of NaN values in matrix

disp(strcat("number of NaN: ",num2str(num_nans)));
disp(strcat("matrix sum: ",num2str(matrix_sum)));
disp(strcat("matrix size: ",num2str(numel(test_matrix))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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
    
clear dconn
newcii.cdata = newdconn;
disp('Saving new Zscored dconn')

save(newcii, [output_cifti_name '.gii'], 'ExternalFileBinary')
disp('Converting Zscored dconn .gii to .nii')
unix([path_wb_c ' -cifti-convert -from-gifti-ext ' output_cifti_name '.gii ' output_cifti_name '.dconn.nii' ]);
disp('Removing .gii')
unix(['rm -f ' output_cifti_name '.gii']);
unix(['rm -f ' output_cifti_name '.dat']);

output_cifti_name = [output_cifti_name '.dconn.nii']; %return the output name with the extension.
end