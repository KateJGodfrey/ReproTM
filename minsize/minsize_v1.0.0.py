import os
import nibabel as nb
import numpy as np
import scipy.io as sio  # for accessing matfiles
from scipy import stats
import argparse

# author: Kate J. Godfrey
# email:  godfreykatej@gmail.com 
# github: https://github.com/KateJGodfrey/ReproTM

# this script identifies network clusters below a minimum size
# and assigns cluster vertices/greyordinates to the mode of neighbors

# get arguements from script call
parser=argparse.ArgumentParser(description='Assign Network Clusters Below Minsize to Mode of Neighbors')
parser.add_argument('--dscalar_infile',required=True,help="path to input dscalar")
parser.add_argument('--dscalar_outfile',required=True,help="path to output dscalar")
parser.add_argument('--minsize',required=True,type=int,help="minimum cluster size")

# get variables from arguements
args = parser.parse_args()
dscalar_infile = args.dscalar_infile
dscalar_outfile = args.dscalar_outfile
minsize = args.minsize

# get code directory
dir_minsize = os.path.dirname(os.path.abspath(__file__))

# get support file directory
dir_support = dir_minsize + '/support_files'

# load data from the dscalar
subject_cii = nb.load(dscalar_infile)
network_assigns = subject_cii.get_fdata()[0]

# get support files
if network_assigns.shape[0] == 59412:
    neighbors_infile = dir_support + '/node_neighbors_59412.mat'
    neighbors = sio.loadmat(neighbors_infile)
    neighbors = np.array(neighbors.get('neighbors'))
    neighbors = neighbors - 1 # subtract one to account for python 0 indexing
elif network_assigns.shape[0] == 91282:
    neighbors_infile = dir_support + '/node_neighbors_91282.mat'
    neighbors = sio.loadmat(neighbors_infile)
    neighbors = np.array(neighbors.get('neighbors'))
    neighbors = neighbors - 1 # subtract one to account for python 0 indexing
else:
    print('something went wrong with input CIFTI of network assignments')
    print('number of network assignments is not 59412 (cortex) or 91282 (whole brain)')
    print('please double check network assignment dscalar')

# new vector to be modified
network_assigns_out = network_assigns

# get total number of unique networks
networks = np.unique(network_assigns)

# perform the operation for each network
for network in networks:
    
    # get vertices assigned to this network
    networkverts = np.where(network_assigns == network)[0]
    
    # get a metric to store unique cluster ids
    clusteredmetric = np.zeros(network_assigns.shape[0])
    
    # step 1: total unique clusters calculation
    for vertex in networkverts:
        
        # get the neighbors of this vertex
        vertexneighbors = neighbors[vertex,:]
        vertexneighbors = vertexneighbors[~np.isnan(vertexneighbors)]   # remove NaN neighbors
        
        # find which neighbors have the same network assignment
        vertexneighbors_thisnetwork = [np.intersect1d(networkverts,vertexneighbors)][0]
        
        # find whether neighbors already have a cluster number
        uniqueneighborvals = np.unique(clusteredmetric[vertexneighbors_thisnetwork])
        uniqueneighborvals = uniqueneighborvals[uniqueneighborvals != 0]
        uniqueneighborvals = uniqueneighborvals.astype(int) 
        
        # give each cluster a unique identifier
        # if no cluster numbers, assign to vertex number
        if len(uniqueneighborvals) == 0:
            cluster = vertex + 1
            clusteredmetric[vertexneighbors_thisnetwork] = cluster
        # if only one cluster number, assign to unique value
        elif len(uniqueneighborvals) == 1:
            cluster = uniqueneighborvals[0]
            clusteredmetric[vertexneighbors_thisnetwork] = cluster
        # if multiple cluster numbers, merge across unique values
        else:
            for valuenum in range(1, len(uniqueneighborvals)):
                clusteredmetric[clusteredmetric == uniqueneighborvals[valuenum]] = uniqueneighborvals[0]
    
    # get the total number of unique clusters
    uniqueclustervals = np.unique(clusteredmetric)
    uniqueclustervals = uniqueclustervals[uniqueclustervals != 0]
    
    # step 2: check if clusters are smaller than minsize and apply mode of neighbors
    for clusternum in uniqueclustervals:
        
        # check how many vertices belong in the cluster
        verts_in_cluster = np.count_nonzero(clusteredmetric == clusternum)
        
        # if this cluster is below the minimum size, 
        if verts_in_cluster < minsize:
            # get the neighbors of this cluster
            neighborverts = neighbors[np.where(clusteredmetric == clusternum),1:][0]
            # get list of neighboring vertices that are unique 
            neighborverts = np.unique(neighborverts)
            # get list of neighbors that are part of the cluster border
            borderverts = np.setdiff1d(neighborverts,np.where(clusteredmetric == clusternum))
            borderverts = borderverts[~np.isnan(borderverts)]
            if network_assigns.shape[0] == 59412:
                borderassigns = network_assigns_out[borderverts[borderverts < 59412]]; 
            # get the assignment of the border verts
            if network_assigns.shape[0] == 91282:
                borderassigns = network_assigns_out[borderverts]; 
            # calculate the mode network assignment of these bordering vertices
            mode_neighborval = stats.mode(borderassigns)[0][0]
            # assign vertices in this cluster to the mode of their neighbors
            network_assigns_out[np.where(clusteredmetric == clusternum)] = mode_neighborval

# save dscalar
dscalar_cii = nb.load(dscalar_infile)
scalar_axis = nb.cifti2.ScalarAxis([]) 
new_img = nb.Cifti2Image(network_assigns_out.reshape(1,network_assigns_out.shape[0])
                        ,header=dscalar_cii.header
                        ,nifti_header=dscalar_cii.nifti_header
                        ,extra=scalar_axis)
new_img.to_filename(dscalar_outfile)
del dscalar_cii, new_img