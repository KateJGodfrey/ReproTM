#!/bin/bash
module load workbench/1.5.0

subject_id=0A4P0LWM
data_dir=/projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/data
session_id=00A
workbench_command=`which wb_command`

INPUT_HDF5=${data_dir}/derivatives/xcp_d_v0.12.0_unstable//sub-${subject_id}/ses-${session_id}/func/sub-${subject_id}_ses-${session_id}_task-rest_desc-abcc_qc.hdf5
OUTPUT_MAT=${data_dir}/derivatives/xcp_d_v0.12.0_unstable//sub-${subject_id}/ses-${session_id}/func/sub-${subject_id}_ses-${session_id}_task-rest_desc-abcc_qc.mat

if [ ! -f "$INPUT_HDF5" ]; then
    echo "Error: Input file '$INPUT_HDF5' not found!"
    exit 1
fi

echo "Converting HDF5 to MAT format..."
python3 /projects/standard/midb_abcd/shared/ABCC/code/precision_mapping_via_template_matching/convert_hdf5_to_matlab.py "$INPUT_HDF5" "$OUTPUT_MAT"
echo "Conversion complete!"

# set up preliminaries for cifti_connectivity
minutes=none
TR=0.8
dtseries_filename=${data_dir}/derivatives/xcp_d_v0.12.0_unstable/sub-${subject_id}/ses-${session_id}/func/sub-${subject_id}_ses-${session_id}_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii
# motion_filename=${data_dir}/derivatives/xcp_d_v0.12.0_unstable//sub-${subject_id}/ses-${session_id}/func/sub-${subject_id}_ses-${session_id}_task-rest_desc-filtered_motion_mask.mat
motion_filename=${OUTPUT_MAT}
cifti_pconn_output_directory=${data_dir}/derivatives/xcp_d_v0.12.0_unstable/cifti_connectivity_outputs/sub-${subject_id}/ses-${session_id}
cifti_conn_wrapper=/projects/standard/midb-ig/shared/projects/ABCD/cifti_conn/cifti-connectivity/cifti_conn_wrapper.py
matlab_compiler=/projects/standard/midb-ig/shared/projects/ABCD/cifti_conn/MCR_V91/v91


# Create cifti pconn output folder if it does not exist
if [ -d ${cifti_pconn_output_directory} ]; then
    rm -rf ${cifti_pconn_output_directory}
    mkdir -p ${cifti_pconn_output_directory}
else
    mkdir -p ${cifti_pconn_output_directory}
fi

# run cifti_connectivity for generating dconns
${cifti_conn_wrapper} --keep-conn-matrices -mre ${matlab_compiler} -wb ${workbench_command} -outliers --minutes 10 --motion ${motion_filename} ${dtseries_filename} ${TR} ${cifti_pconn_output_directory} matrix

# run cifti_connectivity for generating dconns with 10 minutes of data without motion regression
# ${cifti_conn_wrapper} --keep-conn-matrices -mre ${matlab_compiler} -wb ${workbench_command} -outliers --minutes 10 ${dtseries_filename} ${TR} ${cifti_pconn_output_directory} matrix