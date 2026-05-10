#!/bin/bash -l

#SBATCH -J ReproTM
#SBATCH --ntasks=2
#SBATCH --tmp=120gb
#SBATCH --mem=120gb
#SBATCH -t 2:00:00
#SBATCH -p your_partition
#SBATCH -o ReproTM_%j.out
#SBATCH -e ReproTM_%j.err

SUB=${1}
SES=${2}
TASK=${3}

# set directories and paths
# what directory are you in?
dir_in=/path/to/your/directory
# what directory are your input dconns in?
dir_data=/path/to/your/input/dconns
# where do you want to save your template matching outputs?
dir_out=/path/to/your/output/folder

# update output log with script details
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "sub: sub-${SUB}"
echo "ses: ses-${SES}"
echo "task: task-${TASK}"
echo ""

# make sure output folder exists, if not create it
if [ ! -d ${dir_out}/sub-${SUB}/ses-${SES} ]; then
    mkdir -p ${dir_out}/sub-${SUB}/ses-${SES}/func
fi

# input filepaths
# dconn_infile: CIFTI input dconn
dconn_infile=${dir_data}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_stat-zscored.dconn.nii
# template_infile: full path to template, should match the units of input dconn
template_infile=/*/ReproTM/ReproTM/support_files/network_seedmap_templates/choice_of_template.mat
# template_networks: names of your networks, should match template_infile
template_networks="DMN Vis FP NaN DAN NaN VAN Sal AMN SMd SMl Aud Tpole MTL PMN PON NaN SCAN"

# output filepaths
# mat_outfile: ReproTM intermediates including vertex-to-network similarity values and winner-take-all network assignments
mat_outfile=${dir_out}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_ReproTM.mat
# dscalar_outfile: CIFTI output scalar with winner-take-all network assignments
# dsclarSCANrefined_outfile: CIFTI output scalar with winner-take-all network assignments following SCAN network refinement
dscalar_outfile=${dir_out}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_ReproTM.dscalar.nii
dscalarSCANrefined_outfile=${dir_out}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_ReproTM_refineSCAN.dscalar.nii

# user specified options
# --surface_only: optional, run template matching for cortex only.
# --template_thresholding: optional, threshold template seedmap to top connections above threshold.
# --template_minthreshold: if --template_thresholding, what do you want to set as the minimum threshold?
# --refineSCAN: optional, re-run template matching with a higher minimum threshold for SCAN, SMl and SMd networks.
# --refineSCAN_minthreshold: if --refineSCAN, what threshold will be applied to SCAN, SMl and SMd networks?

# path to ReproTM
code_path=/*/ReproTM/ReproTM/ReproTM_v1.0.0.py

# run python
python -u ${code_path} --dconn_infile ${dconn_infile} --template_infile ${template_infile} --template_networks ${template_networks} --template_thresholding --template_minthreshold 1 --refineSCAN --refineSCAN_minthreshold 3 --mat_outfile ${mat_outfile} --dscalar_outfile ${dscalar_outfile} --dscalarSCANrefined_outfile ${dscalarSCANrefined_outfile}

echo "script finished at at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"
