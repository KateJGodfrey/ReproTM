#!/usr/bin/env python3
"""
Append dscalar data to .mat file under 'new_subject_labels_refined_cleaned'
"""

import argparse
import scipy.io as sio
import numpy as np
import nibabel as nib
import sys


def append_dscalar_to_mat(dscalar_filepath, mat_filepath, output_filepath=None):
    """
    Load dscalar data from a .dscalar.nii file and append to a .mat file 
    under 'new_subject_labels_refined_cleaned'
    
    Parameters:
    -----------
    dscalar_filepath : str
        Path to the .dscalar.nii file
    mat_filepath : str
        Path to the input .mat file
    output_filepath : str, optional
        Path to save the updated .mat file. If None, overwrites the input file.
    
    Returns:
    --------
    None
    """
    # Load the dscalar file
    print(f"Loading dscalar from: {dscalar_filepath}")
    try:
        dscalar_img = nib.load(dscalar_filepath)
        dscalar_data = dscalar_img.get_fdata()
    except Exception as e:
        print(f"Error loading dscalar file: {e}")
        sys.exit(1)
    
    print(f"Loaded dscalar shape: {dscalar_data.shape}")
    print(f"Original data type: {dscalar_data.dtype}")
    
    # Ensure it's 2D and get the data (usually shape is (1, 91282) or similar)
    if dscalar_data.ndim == 1:
        dscalar_data = dscalar_data.reshape(1, -1)
    elif dscalar_data.ndim == 2 and dscalar_data.shape[0] > 1:
        # If multiple rows, take the first one or handle as needed
        print(f"Warning: dscalar has {dscalar_data.shape[0]} rows, using first row only")
        dscalar_data = dscalar_data[0:1, :]
    
    # Convert to int64
    dscalar_data = dscalar_data.astype(np.int64)
    print(f"Converted to int64")
    
    # Verify dimensions
    if dscalar_data.shape[1] != 91282:
        print(f"Error: Expected 91282 columns, got {dscalar_data.shape[1]}")
        sys.exit(1)
    
    # Load the .mat file (or create empty dict if file doesn't exist)
    try:
        mat_data = sio.loadmat(mat_filepath)
        print(f"Loaded existing .mat file: {mat_filepath}")
    except FileNotFoundError:
        mat_data = {}
        print(f"File not found. Creating new .mat file.")
    except Exception as e:
        print(f"Error loading .mat file: {e}")
        sys.exit(1)
    
    # Check if field already exists
    if 'new_subject_labels_eta_refineSCAN_minsize30' in mat_data:
        existing_data = mat_data['new_subject_labels_eta_refineSCAN_minsize30']
        print(f"Existing data shape: {existing_data.shape}")
        print(f"Existing data type: {existing_data.dtype}")
        
        # Ensure existing data is also int64
        if existing_data.dtype != np.int64:
            existing_data = existing_data.astype(np.int64)
            print(f"Converted existing data to int64")
        
        print(f"Existing data shape: {existing_data.shape} - will be overwritten")
        updated_data = dscalar_data
    else:
        # Create new field
        updated_data = dscalar_data
        print(f"Created new field with shape: {updated_data.shape}")
    
    # Ensure final data is int64
    updated_data = updated_data.astype(np.int64)
    print(f"Final data type: {updated_data.dtype}")
    
    # Update the mat_data
    mat_data['new_subject_labels_eta_refineSCAN_minsize30'] = updated_data
    
    # Determine output path
    if output_filepath is None:
        output_filepath = mat_filepath
    
    # Save the updated .mat file
    try:
        sio.savemat(output_filepath, mat_data)
        print(f"Successfully saved to: {output_filepath}")
    except Exception as e:
        print(f"Error saving .mat file: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Append dscalar data to .mat file under new_subject_labels_refined_cleaned (as int64)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Append dscalar to mat file (overwrites)
  %(prog)s --dscalar seedmap.dscalar.nii --mat data.mat
  
  # Append and save to new file
  %(prog)s --dscalar seedmap.dscalar.nii --mat input.mat --output output.mat
  
  # Short flags
  %(prog)s -d seedmap.dscalar.nii -m data.mat -o output.mat
        """
    )
    
    parser.add_argument(
        '--dscalar', '-d',
        required=True,
        type=str,
        help='Path to the input .dscalar.nii file'
    )
    
    parser.add_argument(
        '--mat', '-m',
        required=True,
        type=str,
        help='Path to the input .mat file (will be created if does not exist)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Path to the output .mat file (optional, defaults to overwriting input mat file)'
    )
    
    args = parser.parse_args()
    
    # Run the function
    append_dscalar_to_mat(
        dscalar_filepath=args.dscalar,
        mat_filepath=args.mat,
        output_filepath=args.output
    )


if __name__ == "__main__":
    main()