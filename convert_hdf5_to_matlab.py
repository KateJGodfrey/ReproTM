#!/usr/bin/env python3
import h5py
import scipy.io as sio
import numpy as np
import sys

def convert_hdf5_to_mat(input_file, output_file):
    # Open your HDF5 file
    with h5py.File(input_file, 'r') as f:
        dcan_motion = f['dcan_motion']
        
        # Get all subject keys
        subject_keys = list(dcan_motion.keys())
        print(f"Number of subjects: {len(subject_keys)}")
        
        # Create a cell array for MATLAB (using object array in NumPy)
        # Changed shape to (1, 101) instead of (101,)
        cell_array = np.empty((1, len(subject_keys)), dtype=object)
        
        # Populate each cell with the subject's data
        for i, subject_key in enumerate(subject_keys):
            subject_group = dcan_motion[subject_key]
            
            # Create a dictionary to hold all datasets for this subject
            subject_data = {
                'subject_id': subject_key,
                'FD_threshold': float(subject_key.split('_')[-1])  # Add FD_threshold with the same value
            }
            
            # Read all datasets in this subject's group
            for dataset_name in subject_group.keys():
                dataset = subject_group[dataset_name]
                
                # Rename binary_mask to frame_removal
                output_name = 'frame_removal' if dataset_name == 'binary_mask' else dataset_name
                
                # Read scalar or array data
                if dataset.shape == ():
                    subject_data[output_name] = dataset[()]
                else:
                    subject_data[output_name] = dataset[:]
            
            # Store in cell array - note the [0, i] indexing
            cell_array[0, i] = subject_data
        
        # Save to .mat file
        sio.savemat(output_file, {'motion_data': cell_array}, oned_as='row')
    
    print(f"\nSuccessfully saved to {output_file}")
    
    # Verify
    mat_data = sio.loadmat(output_file)
    motion_data = mat_data['motion_data']
    print(f"Verified: motion_data shape = {motion_data.shape}")

if __name__ == "__main__":
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'your_file.hdf5'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'dcan_motion_data.mat'
    convert_hdf5_to_mat(input_file, output_file)