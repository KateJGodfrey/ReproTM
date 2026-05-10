function [outname] = clean_dscalars_by_size(dscalarwithassignments,manualset,groupnetworksfile,dostripes,mincol,minsize,orig_parcelsfile,make_consensus,assign_unassigned,remove46)
%'${dscalar_infile}',[],[],[],[],30,[],0,0,1
%consensus_maker_knowncolors(regularized_ciftifile,[manualset],[groupnetworksfile],[dostripes],[mincol],[minsize],[orig_parcelsfile])
% clean_dscalars_by_size([cifti_output_folder '/' output_cifti_name '_overlap_' overlap_method '.dtseries.nii'],[],[],[],[],30,[],0,0,1,0);
%This function cleans up networks.
%Hardcode
%make_consensus = 1;
%dscalarwithassignments='/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/sub-10227/ses-combined/func/sub-10227_ses-combined_task-crossMENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.dscalar.nii'
%minsize=30;
%assign_unassigned=0;

dscalarwithassignments = dscalarwithassignments;
% disp(dscalarwithassignments);
% class(dscalarwithassignments);
if exist(dscalarwithassignments) == 0
    NOTE = ['subject dscalar does not exist']
    disp(dscalarwithassignments);
    return
else
    disp('all series files exist continuing ...')
end

network_assignment_filetype = strsplit(dscalarwithassignments, '.');
cifti_type = char(network_assignment_filetype(end-1));
if strcmp('dtseries',cifti_type) == 1
    overlap =1;
else
    overlap =0;
end

%% Adding paths for this function
this_code = which('clean_dscalars_by_size');
support_folder='/projects/standard/faird/shared/code/internal/analytics/compare_matrices_python/code/clean_dscalars/support_files';
addpath(genpath(support_folder));

% where is your library for reading and writing cifti files?
cifti_codebase='/projects/standard/faird/shared/code/external/utilities/cifti-matlab';
addpath(genpath(cifti_codebase));


if ~exist('mincol','var') || isempty(mincol)
    mincol = 1;
end

if ~exist('minsize','var') || isempty(minsize)
    minsize = 0;
end

% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
all_color_values = [1:100];


%% Start code

regularized_ciftifile = dscalarwithassignments;
cifti_data =ciftiopen(regularized_ciftifile); 
assigns = cifti_data.cdata;
assigns(assigns<0) = 0; %Assignments that are eqaul to -1 (unassigned), set the to 0.
assigns(isnan(assigns)) = 0; % Set nans to 0.

out = assigns;
temp_out = out;
temp_all =  temp_out;

if size(temp_out,1) == 91282
    cifti_template_infile='/projects/standard/faird/shared/code/internal/analytics/compare_matrices_python/code/clean_dscalars/support_files/91282_Greyordinates.dscalar.nii';
elseif size(temp_out,1) == 59412
    cifti_template_infile='/projects/standard/faird/shared/code/internal/analytics/compare_matrices_python/code/clean_dscalars/support_files/91282_Greyordinates_surf_only.dscalar.nii'
else
    disp('Something went wrong, check input cifti size')
end

if logical(minsize)
    % Clean up tiny twinspieces
    if size(temp_out,1) == 91282
        if exist([support_folder '/node_neighbors_91282.mat'],'file') == 2
            load([support_folder '/node_neighbors_91282.mat']);
        else
            path_to_neighbors_file = [support_folder filesep 'node_neighbors_59412.mat'];
            neighbors = cifti_neighbors_dcan(regularized_ciftifile,[],[],path_to_neighbors_file);
        end
    else
        if exist([support_folder '/node_neighbors_59412.mat'],'file') == 2
            load([support_folder '/node_neighbors_59412.mat']);
        else
            path_to_neighbors_file = [support_folder filesep 'node_neighbors.txt'];
            neighbors = cifti_neighbors_dcan(regularized_ciftifile,[],[],path_to_neighbors_file);
        end
    end
    for j = 1:size(out,2)
        %if j==4 || j ==6
        if j==3 | j==7 | j==9 | j==12 | j==14 | j==15 | j==16 | j==17
            continue
        end
        allcolors= unique(out(:,j));
        if assign_unassigned == 1
            %allcolors = [allcolors(2:end)' allcolors(1)']'; % move unsigned networks to the end, so they're likely to be already assinged.
        else %leave unsigned vertices unassinged.
            allcolors(allcolors<=0) = [];
        end
        temp_out = temp_all(:,j);
        for color = allcolors(:)'
            clusteredmetric = zeros(size(temp_out));
            thiscolorverts = find(temp_out==color);
            for vertex = thiscolorverts'
                %find the neighbors of this vertex
                vertexneighbors = neighbors(vertex,:);
                vertexneighbors(isnan(vertexneighbors)) = [];
                %find which of those neighbors also pass the thresholds
                vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
                %find if those neighbors have already been assigned different cluster values
                uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
                uniqueneighborvals(uniqueneighborvals==0) = [];
                %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
                if isempty(uniqueneighborvals)
                    clusteredmetric(vertexneighbors_thiscolor) = vertex;
                    %if there is only one previous cluster identifier present, make all the neighbors that value
                elseif length(uniqueneighborvals)==1
                    clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
                    %if there are multiple cluster identifier values in the neighborhood, merge them into one
                else
                    for valuenum = 2:length(uniqueneighborvals)
                        clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
                        if color == 0
                            if assign_unassigned ==1 % added by Robert H.
                                clusteredmetric(vertex) = uniqueneighborvals(1); % don't forget to assign the vertex the same assingment.
                            else
                            end
                        else
                        end
                    end
                end
            end
            uniqueclustervals = unique(clusteredmetric);
            uniqueclustervals(uniqueclustervals==0) = [];
            for clusternum = uniqueclustervals' %added by Robert
                %                disp(['cluster ID ' num2str(clusternum)])
                verts_in_cluster = nnz(clusteredmetric==clusternum);
                %                     disp([num2str(verts_in_cluster) ' vertices in cluster'])
                %                end
                %if nnz(clusteredmetric==clusternum) < minsize
                if verts_in_cluster < minsize
                    neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
                    neighborverts(isnan(neighborverts)) = [];
                    borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
                    if size(temp_out,1) ~= 91282
                        large_val=find(borderverts>59412);
                        borderverts(large_val)=NaN;
                        if size(large_val)==size(borderverts) %added by Robert 05/31
                            disp('For some reason all the indices of the neighbors for this cluster are larger than the size of this dscalar');
                            disp(['clusternum is: ' num2str(clusternum) '. Using neighbors neighbors...']);
                            extended_neighborverts = unique(neighbors(clusternum-1,2:end)); %grab a neighbors of the grayordinate that is 1 value lower.  Not a great fix admittedly -RH
                            extended_neighborverts(isnan(extended_neighborverts)) = [];
                            borderverts = setdiff(extended_neighborverts,find(clusteredmetric==clusternum));
                            large_val=find(borderverts>59412);
                            borderverts(large_val)=NaN;
                        end
                        borderverts=borderverts(~isnan(borderverts));
                    end
                    %borderverts(temp_out(borderverts)<1) = [];
                    %added by robert
                    if assign_unassigned ==1
                        %borderverts(temp_out(borderverts)) = [];
                        bordererassigns = temp_out(borderverts);
                        int_clusterassings = bordererassigns > 0;
                        mode_neighborval = mode(bordererassigns(int_clusterassings));
                        if mode_neighborval ==0
                            disp(clusternum)
                        end
                    else
                        if overlap ==0 % allow networks to get an assingment of 0, but only for overlapping networks.
                            borderverts(temp_out(borderverts)<1) = [];
                            %mode_neighborval = mode(temp_out(isassinged(borderverts)));
                        else
                            borderverts(temp_out(borderverts)<0) = [];
                        end
                        mode_neighborval = mode(temp_out(borderverts));
                        %                             if isnan(mode_neighborval) ==1
                        %                                 mode_neighborval =0;
                        %                             end
                    end
                    %Grab the next value
                    %mode_neighborval = mode(temp_out(borderverts));
                    temp_out(clusteredmetric==clusternum) = mode_neighborval;
                else % don't forget to assign large clusters to networks (network ==0) to values is still unassigned.
                    if assign_unassigned ==1
                        if color == 0
                            neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
                            neighborverts(isnan(neighborverts)) = [];
                            borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
                            if size(temp_out,1) ~= 91282
                                large_val=find(borderverts>59412);
                                borderverts(large_val)=NaN;
                                if size(large_val)==size(borderverts) %added by Robert 05/31
                                    disp('For some reason all the indices of the neighbors for this cluster are larger than the size of this dscalar');
                                    disp(['clusternum is: ' num2str(clusternum) '. Using neighbors neighbors...']);
                                    extended_neighborverts = unique(neighbors(clusternum-1,2:end)); %grab a neighbors of the grayordinate that is 1 value lower.  Not a great fix admittedly -RH
                                    extended_neighborverts(isnan(extended_neighborverts)) = [];
                                    borderverts = setdiff(extended_neighborverts,find(clusteredmetric==clusternum));
                                    large_val=find(borderverts>59412);
                                    borderverts(large_val)=NaN;
                                end
                                borderverts=borderverts(~isnan(borderverts));
                            end
                            bordererassigns = temp_out(borderverts);
                            int_clusterassings = bordererassigns > 0;
                            mode_neighborval = mode(bordererassigns(int_clusterassings));
                            temp_out(clusteredmetric==clusternum) = mode_neighborval;
                        else
                        end
                    else
                    end
                end
            end
        end
        all_temp_out(:,j) = temp_out;
    end
    out = temp_out;
end
if overlap == 1
    cifti_data.data = all_temp_out;
else
    cifti_data.cdata = out;
end

% if overlap == 0
%     if ~exist('cifti_data.mapname')
%         cifti_data.mapname = {'Column number ' num2str(mincol)};
%         cifti_data.dimord = 'scalar_pos';
%     else
%         cifti_data.mapname = cifti_data.mapname(mincol);
%     end
% else
%     if ~exist('cifti_data.mapname')
%         cifti_data.mapname = {'Column number ' num2str(mincol)};
%         cifti_data.dimord = 'pos_time';
%     else
%         cifti_data.mapname = cifti_data.mapname(mincol);
%     end
% end

dotsloc = strfind(regularized_ciftifile,'.');
basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
outname = [ basename '_recolored' '_minsize' num2str(minsize) '.dscalar.nii']; 

cifti_template = cifti_read(cifti_template_infile);
cifti_template.cdata = out;

cifti_write(cifti_template,outname);

% if overlap ==0
%     set_cifti_powercolors([outname '.dscalar.nii'])
% else
%     set_cifti_powercolors([outname '.dtseries.nii'])
% end
% if exist('orig_parcelsfile') && ~isempty(orig_parcelsfile)
%     parcels = ft_read_cifti_mod(orig_parcelsfile);
%     parcels = parcels.data;
%     parcels((length(out)+1):end) = [];
%     IDs = unique(parcels); IDs(IDs<1) = [];
%     outtext = zeros(length(IDs),1);
%     outtext_bycol = zeros(length(IDs),size(all_recolored,2));
%     for IDnum = 1:length(IDs)
%         outtext(IDnum) = mode(out(parcels==IDs(IDnum)));
%         outtext_bycol(IDnum,:) = mean(all_recolored(parcels==IDs(IDnum),:),1);
%     end
%     dlmwrite([outname '.txt'],outtext,'delimiter',' ')
%     dlmwrite([basename '_allcolumns_recolored.txt'],outtext_bycol,'delimiter',' ')
% end

disp('Done running clean dscalars code.')