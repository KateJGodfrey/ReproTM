#!/bin/bash -l

#SBATCH -J minsize
#SBATCH -c 1
#SBATCH --tmp=5gb
#SBATCH --mem=5gb
#SBATCH -t 5:00
#SBATCH -p your_partition
#SBATCH -o minsize_%j.out
#SBATCH -e minsize_%j.err

# author: Kate J. Godfrey 
# github: https://github.com/KateJGodfrey/ReproTM

# this script can be used to create input variables and call minsize
# minsize identifies network clusters below a minimum size (in greyordinates)
# and assigns cluster vertices/greyordinates to the mode network assignment of neighbors

# user specified terminal inputs
SUB=${1}
SES=${2}
TASK=${3}

# set paths for input and output dscalars
dir_data=/path/to/your/input/dscalars
dscalar_infile=${dir_data}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_ReproTM.dscalar.nii
dscalar_outfile=${dir_data}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}_ReproTM_minsize30.dscalar.nii
minsize=30

# python code 
code_path=/*/code/minsize/minsize_v1.0.0.py

# run script
python ${code_path} --dscalar_infile ${dscalar_infile} --dscalar_outfile ${dscalar_outfile} --minsize 30
