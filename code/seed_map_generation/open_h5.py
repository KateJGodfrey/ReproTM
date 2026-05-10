import h5py

# Check if the file is valid HDF5
try:
    with h5py.File('/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/compare_matrices_python/code/seed_map_generation/all_data/python_seed_map_refractor/seedmaps_all_data_dtsereies.mat', 'r') as f:
        print("File is valid HDF5")
        print("Variables:", list(f.keys()))
except:
    print("File is NOT valid HDF5")