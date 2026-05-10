import os
import sys
import nibabel as nb
import scipy.io as sio       # scipy is for accessing mat files
import numpy as np
import time
import argparse

# get the current time
totaltimer = time.time()

# get code directory
dir_ReproTM = os.path.dirname(os.path.abspath(__file__))

# get support file directory
dir_support = dir_ReproTM + '/support_files'

# get args
parser=argparse.ArgumentParser(description='Reproducible Template Matching')
parser.add_argument('--dconn_infile',required=True,help="path to input dconn")
parser.add_argument('--template_infile',required=True,help="path to input template matfile")
parser.add_argument('--template_networks',required=True,nargs="*",help="network names, should match input template")
parser.add_argument('--surface_only',required=False,action='store_true'
                    ,help="perform template matching for cortex only?")
parser.add_argument('--template_thresholding',required=False,action='store_true',
                    help="threshold input template?")
parser.add_argument('--template_minthreshold',required=False,type=int,help="if --template_thresholding, set minimum threshold")
parser.add_argument('--mat_outfile',required=True,help="path for output matfile")
parser.add_argument('--dscalar_outfile',required=True,help="path for output dscalar with network assignments")
parser.add_argument('--refineSCAN',required=False,action='store_true'
                    ,help="optional: re-run template matching to refine SCAN network")
parser.add_argument('--refineSCAN_minthreshold',required=False,type=int,help="if --refineSCAN, set minimum threshold")
parser.add_argument('--dscalarSCANrefined_outfile',required=False,help="if --refineSCAN, path for SCAN refined output dscalar")
args = parser.parse_args()

# get vars
# input dconn
dconn_infile = args.dconn_infile
# input template
template_infile = args.template_infile
# input template network names
template_networks = args.template_networks
# output matlab file
mat_outfile = args.mat_outfile
# output dscalar files
dscalar_outfile = args.dscalar_outfile
# surface only or whole brain
surface_only=args.surface_only
# thresholding options
template_thresholding=args.template_thresholding
if args.template_thresholding and args.template_minthreshold is None:
    parser.error("--template_minthreshold is required for --template_thresholding")
    sys.exit()
if args.template_thresholding:
    template_minthreshold=args.template_minthreshold
# refinement options
refineSCAN=args.refineSCAN
if refineSCAN:
    # if refineSCAN, check for additional arguments
    if args.refineSCAN_minthreshold is None:
        parser.error("--refineSCAN_minthreshold is required for --refineSCAN")
        sys.exit()
    if args.dscalarSCANrefined_outfile is None:
        parser.error("--dscalarSCANrefined_outfile is required for --refineSCAN")
        sys.exit()
    # if refineSCAN, assign variables
    refineSCAN_networks = ['SMd','SMl','SCAN']
    refineSCAN_minthreshold=args.refineSCAN_minthreshold
    dscalarSCANrefined_outfile=args.dscalarSCANrefined_outfile


# can hardcode paths to necessary directories, files and code
# this makes it easier to debug the code and run on command line
hardcode=False
if hardcode:
    # workbench path
    wb_command='/path/workbench'
    # input file paths
    dir_in='/path/working_dir'
    dconn_infile = '/path/input_dconn_stat-zscored.dconn.nii'
    template_infile='/*/ReproTM/support_files/network_seedmap_templates/input_template.mat'
    # output file paths
    mat_outfile='/path/output_ReproTM.mat'
    dscalar_outfile='/path/output_ReproTM.dscalar.nii'
    dscalarSCANrefined_outfile='/path/output_ReproTM_refineSCAN.dscalar.nii'
    # user specified options
    surface_only=True
    template_thresholding=True
    template_minthreshold=float(1.00)
    refineSCAN=False
    refineSCAN_minthreshold = float(3.00)
    template_networks = ['DMN','Vis','FP','NaN','DAN','NaN','VAN','Sal','AMN','SMd','SMl','Aud','Tpole','MTL','PMN','PON','NaN','SCAN']

# if thresholding, make sure thresholds are numeric
if template_thresholding:
    template_minthreshold = float(template_minthreshold)

if refineSCAN:
    refineSCAN_minthreshold = float(refineSCAN_minthreshold)

# determine template to use for output dscalar based on surface_only
if surface_only:
    dscalar_template = dir_support + '/cifti_templates/59412_SurfaceVerticesCifti2.dscalar.nii'

if not surface_only:
    dscalar_template = dir_support + '/cifti_templates/91282_GreyordinatesCifti2.dscalar.nii'

# check that necessary files exist
file_checker=True
if file_checker:
    # check for input dconn
    if os.path.isfile(dconn_infile) == False:
        print('your dconn path doesnt exist')
        print('\tinput dconn file path provided: ' + dconn_infile)
        sys.exit('')
    # check for input network template seedmaps
    if os.path.isfile(template_infile) == False:
        print('your template path doesnt exist')
        print('\tinput template file path provided: ' + template_infile)
        sys.exit('')
    # check for dscalar template in support folder
    if os.path.isfile(dscalar_template) == False:
        print('your dscalar template path doesnt exist')
        print('\tReproTM needs cifti template for network assignment .dscalar.nii')
        print('\tinput dscalar template path: ' + dscalar_template)
        print('\tcheck whether your codebase has the necessary support folder')

load_data=True
if load_data:
    print('loading data')
    print('\t')
    steptimer = time.time()
    
    print('\tloading the correlation matrix')
    subject_cii = nb.load(dconn_infile)
    corr_mat_full = subject_cii.get_fdata()  # get dconn out of input cifti
    greyordinates = corr_mat_full.shape[0]  # get total number of greyordinates from input cifti
    del subject_cii # save memory
    print('\tnumber of greyordinates: ' + str(greyordinates))
    print('\t ')
    
    print('\tloading the template')
    cifti_template = sio.loadmat(template_infile)
    cifti_template_mat_full = cifti_template.get("seed_matrix")
    print('\tnumber of template greyordinates: ' + str(cifti_template_mat_full.shape[0]))
    print('\t')
    
    if surface_only:
        print('\tremoving subcortical data')
        corr_mat_full = corr_mat_full[0:59412,0:59412] # 59412 surface greyordinates
        greyordinates = corr_mat_full.shape[0]
        cifti_template_mat_full = cifti_template_mat_full[0:59412,] 
        print('\tnumber of greyordinates: ' + str(greyordinates))
        print('\tnumber of template greyordinates: ' + str(cifti_template_mat_full.shape[0]))
        print('\t')
    
    if template_thresholding:
        cifti_template_mat_full[cifti_template_mat_full <= template_minthreshold] = 'NaN'
        print('\ttemplate minimum threshold set at: ' + str(template_minthreshold))
        print('\ttemplate greyordinates below minimum set to NaN')
        print('\t')
    if template_thresholding==False:
        print('\tno minimum threshold applied to template')
        print('\t')
        
# check dconn and template maximum values
# will throw an error if template is transformed but dconn is not
if np.max(corr_mat_full) == 1 and np.nanmax(cifti_template_mat_full) > 1:
    print('\t WARNING: INPUT DCONN AND INPUT TEMPLATE SEEM TO DIFFER IN RANGE')
    print('\t if dconn and template are not similarly transformed it will cause errors with template matching')
    print('\t check dconn and template transformations (both r-correlation or z-transformed r-correlation)')
    print('\t script for dconn transformation to z-scored dconn is avilable in the codebase directory')
    sys.exit(1)

# check dconn and template maximum values
# will throw an error if dconn is transformed but template is not
if np.max(corr_mat_full) > 1 and np.max(cifti_template_mat_full) == 1:
    print('\t WARNING: INPUT DCONN AND INPUT TEMPLATE SEEM TO DIFFER IN RANGE')
    print('\t if dconn and template are not similarly transformed it will cause errors with template matching')
    print('\t check dconn and template transformations (both r-correlation or z-transformed r-correlation)')
    sys.exit(1)

# check dconn and template number of greyordinates
# will throw an error if template and dconn differ in number of greyordinates
if (corr_mat_full.shape[0] != cifti_template_mat_full.shape[0]):
    print('\t WARNING: INPUT DCONN AND INPUT TEMPLATE DIFFER IN NUMBER OF GREYORDINATES')
    print('\t if dconn and template have an unequal number of greyordinates will cause errors with template matching')
    print('\t check resolution of inputs (cifti 32k or 164k) or inclusion/exclusion of subcortical greyordinates')
    sys.exit(1)

print('\tloading data complete')
print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
print('\t ')

run_template_matching=True
if run_template_matching:
    
    print('template matching has begun')
    print('\t')
    steptimer = time.time()
    
    # create empty array to store template matched outputs (1 row per greyordinate and 1 column per network)
    # eta_to_template_vox = calculate association to template using sum of squares variance decomposition
    # r_to_template_vox = calculate association to template with linear bivariate r correlation
    eta_to_template_vox = np.zeros((int(len(corr_mat_full)),int(len(template_networks))))
    r_to_template_vox = np.zeros((int(len(corr_mat_full)),int(len(template_networks))))
    
    print('\tcalculating similarity of greyordinates to template')
    print('\t')
    
    # i will loop through each vertex from the correlation matrix
    for i in range(greyordinates):
        
        if i % 5000 == 0: # print a progress update every 5000 voxels
            print('\tcalculating greyordinate: ' + str(i))
            print('\ttemplate matching has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
        for j in range(len(template_networks)): # j will loop through each of the candidate networks to match to
            if template_networks[j] == 'NaN':
                pass
            else:
                # identify greyordinates which are not 'NaN' in cmap and tmap
                # 'NaN' greyordinates specified earlier w/ user-specified tmap threshold
                goodvox = np.where(~np.isnan(corr_mat_full[i,:]) & ~np.isnan(cifti_template_mat_full[:,j]))
                # cmap is a vector, representing is a subset of correlation_matrix for greyordinate i that aren't NaN
                cmap = corr_mat_full[i,goodvox]
                # tmap is vector, representing the template map for network 'j' greyordinates that aren't NaN
                tmap = cifti_template_mat_full[goodvox,j]
                
                # calculation of eta values using within and between sum of squares ##########################################
                # calculate the grand mean across all elements of tmap and cmap
                Mgrand  = ((np.sum(tmap)/(tmap.shape[1])) + (np.sum(cmap)/(cmap.shape[1]))) / 2
                # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                Mwithin = np.mean([tmap, cmap], axis=0)
                # get a sum of squares within
                SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                # get a sum of squares total
                SStot    = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                # eta is calculated as 1 - this proportional variance explained
                eta = (1 - (SSwithin/SStot))
                eta_to_template_vox[i,j] = (1 - (SSwithin/SStot))
                
                # calculation of r values using correlation ######################################################################
                r_to_template_vox[i,j] = np.corrcoef(cmap, tmap)[0,1]
                # delete outputs before the next loop
                del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox
    
    # winner-take-all assignment of greyordinates to networks
    print('\t')
    print('\twinner-take-all assignment of greyordinates to networks')
    print('\t')
    # get new network assignments by using idxmax to return which column has the highest eta
    new_subject_labels = np.nanargmax(eta_to_template_vox,axis=1)
    # get new network assignments by using idxmax to return which column has the highest r 
    new_subject_labels_r = np.nanargmax(r_to_template_vox,axis=1)
    # wrap up this step
    print('\ttemplate matching complete')
    print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
    print('\t')

# optional template matching refinement for SMd, SMl, and SCAN networks with new threshold
if refineSCAN:
    
    print(' ')
    print('template matching optional SCAN refinement has begun')
    print('')
    refinetimer = time.time()
    
    # copy eta results into a new array
    eta_to_template_vox_modified = eta_to_template_vox.copy()
    r_to_template_vox_modified = r_to_template_vox.copy()
    
    # copy winner-take-all results into new array
    new_subject_labels_refined=new_subject_labels.copy()
    new_subject_labels_r_refined=new_subject_labels_r.copy()
    
    # get the seed matrix from before but apply the new threshold
    #cifti_template = sio.loadmat(template_infile)
    #cifti_template = cifti_template.get("seed_matrix")
    cifti_template_mat_full_refined = cifti_template_mat_full.copy()
    cifti_template_mat_full_refined[cifti_template_mat_full_refined <= refineSCAN_minthreshold] = 'NaN'
    
    # refine motor greyordinates
    refine_motor_grays=True
    if refine_motor_grays:
        
        print('\trefinement of motor grays:')
        steptimer=time.time()
        
        # find greyordinates assigned to SCAN refinement networks (SMd, SMl, SCAN) via indexing
        motor_network_idx = [template_networks.index(motor_network) for motor_network in refineSCAN_networks]
        motor_grays = np.transpose(np.where(np.isin(new_subject_labels_refined, motor_network_idx)))
        
        # run template matching looping over all motor_grays and all networks
        for i in range(motor_grays.shape[0]):
            
            if i % 5000 == 0:
                print('\t\tmotor grey refinement has run for: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
                
            for j in range(eta_to_template_vox_modified.shape[1]):
                
                # skip networks that don't exist
                if template_networks[j] == 'NaN':
                    pass
                
                else:
                    
                    # identify voxels which are not 'NaN' in cmap and tmap
                    motor_gray = motor_grays[i].item()
                    goodvox = np.where(~np.isnan(corr_mat_full[motor_gray,:]) & ~np.isnan(cifti_template_mat_full_refined[:,j]))
                    # cmap is a vector, correlation_matrix for greyordinate i that aren't NaN
                    cmap = corr_mat_full[motor_gray,(goodvox)]
                    # tmap is vector, template map for network 'j' greyordinates that aren't NaN
                    tmap = cifti_template_mat_full_refined[goodvox,j]
                    
                    # check if cmap or tmap are zero, and if so return 'nan'
                    if (np.sum(cmap) == 0) and (np.sum(tmap) == 0):
                        eta_to_template_vox_modified[motor_gray,j] = float("nan")
                        r_to_template_vox_modified[i,j] = float("nan")
                    
                    # if cmap and tmap are not zero calculate effect size as usual
                    else:
                        
                        # calculation of eta values using within and between sum of squares ##########################################
                        # calculate the grand mean across all elements of tmap and cmap
                        Mgrand = (np.mean(tmap) + np.mean(cmap)) / 2
                        # calculate the element-wise mean of cmap and tmap to produce a vector of the same length as cmap and tmap
                        Mwithin = np.mean([tmap, cmap], axis=0)
                        # get a sum of squares within
                        SSwithin = (np.sum((tmap-Mwithin)**2)) + (np.sum((cmap-Mwithin)**2))
                        # get a sum of squares total
                        SStot    = (np.sum((tmap-Mgrand)**2)) + (np.sum((cmap-Mgrand)**2))
                        # eta is calculated as 1 - this proportional variance explained
                        eta = (1 - (SSwithin/SStot))
                        eta_to_template_vox_modified[motor_gray,j] = (1 - (SSwithin/SStot))
                        
                        # calculation of r values using correlation ######################################################################
                        r_to_template_vox_modified[motor_gray,j] = np.ma.corrcoef(cmap, tmap)[0,1]
                        # delete outputs before the next loop
                        del cmap, tmap, Mgrand, Mwithin, SSwithin, SStot, goodvox
        
        print('\t\trefinement of motor grays complete')
        print('\t\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
        print('\t\t')
    
    print('\t')
    print('\twinner-take-all re-assignment of greyordinates to networks')
    print('\t')
    new_subject_labels_refined = np.nanargmax(eta_to_template_vox_modified,axis=1)
    new_subject_labels_r_refined = np.nanargmax(r_to_template_vox_modified,axis=1)
    
    print('\ttemplate matching refinement network complete')
    print('\tindividual step took: ' + str(round((time.time()-refinetimer)/60,3)) + ' minutes')
    print('\t')

save_outputs=True
if save_outputs:
    
    print('saving outputs')
    print('\t')
    steptimer=time.time()
    
    print('\tsaving intermediates to .mat file')
    
    # create different dictionaries depent on running refineSCAN
    if refineSCAN:
        
        # b/c python is base zero, add 1 to each element of new_subject_labels to change network labels from 0:17 -> 1:18
        new_subject_labels = np.array([x + 1 for x in new_subject_labels])
        new_subject_labels_refined = np.array([x + 1 for x in new_subject_labels_refined])
        new_subject_labels_r = np.array([x + 1 for x in new_subject_labels_r])
        new_subject_labels_r_refined = np.array([x + 1 for x in new_subject_labels_r_refined])
        
        # combine outputs into a dictionary
        mat_dict={"network_names": np.array(template_networks,dtype=object)
                ,"eta_to_template_vox": eta_to_template_vox
                ,"eta_to_template_vox_modified": eta_to_template_vox_modified
                ,"new_subject_labels_eta": new_subject_labels
                ,"new_subject_labels_eta_refineSCAN": new_subject_labels_refined
                ,"r_to_template_vox": r_to_template_vox
                ,"r_to_template_vox_modified": r_to_template_vox_modified
                ,"new_subject_labels_r": new_subject_labels_r
                ,"new_subject_labels_r_refineSCAN": new_subject_labels_r_refined
                }
        
        # the actual saving of the dictionary with scipy
        sio.savemat(mat_outfile,mat_dict,appendmat=False)
        
    else:
        new_subject_labels = np.array([x + 1 for x in new_subject_labels]).reshape(1,greyordinates)
        new_subject_labels_r = np.array([x + 1 for x in new_subject_labels_r]).reshape(1,greyordinates)
        mat_dict={"network_names": np.array(template_networks,dtype=object)
                ,"eta_to_template_vox": eta_to_template_vox
                ,"new_subject_labels_eta":new_subject_labels
                ,"r_to_template_vox":r_to_template_vox
                ,"new_subject_labels_r":new_subject_labels_r}
        # the actual saving of the dictionary with scipy
        sio.savemat(mat_outfile,mat_dict,appendmat=False)
    
    print('\tsaving intermediates to .mat file complete')
    print('\t')
    
    print('\tsaving subject dscalar file')
    
    # save the template matching winner-take-all networks to CIFTI dscalar
    dscalar_cii = nb.load(dscalar_template)
    scalar_axis = nb.cifti2.ScalarAxis([]) 
    new_img = nb.Cifti2Image(new_subject_labels.reshape(1,greyordinates)
                            ,header=dscalar_cii.header
                            ,nifti_header=dscalar_cii.nifti_header
                            ,extra=scalar_axis)
    new_img.to_filename(dscalar_outfile)
    
    # save the template matching winner-take-all networks following SCAN refinement to CIFTI dscalar, if applicable
    if refineSCAN:
        dscalar_cii = nb.load(dscalar_template)
        scalar_axis = nb.cifti2.ScalarAxis([])
        new_img = nb.Cifti2Image(new_subject_labels_refined.reshape(1,greyordinates)
                                ,header=dscalar_cii.header
                                ,nifti_header=dscalar_cii.nifti_header
                                ,extra=scalar_axis)
        new_img.to_filename(dscalarSCANrefined_outfile)
    
    print('\tsaving subject dscalar file complete')
    print('\t')
    
    # wrap up this step
    print('\tsaving outputs complete')
    print('\tindividual step took: ' + str(round((time.time()-steptimer)/60,3)) + ' minutes')
    print('\t')
