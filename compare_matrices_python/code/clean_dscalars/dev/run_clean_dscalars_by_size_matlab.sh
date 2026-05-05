#!/bin/bash -l

#SBATCH -J clean_dscalar
#SBATCH --ntasks=4
#SBATCH --tmp=5gb
#SBATCH --mem=5gb
#SBATCH -t 0:20:00
#SBATCH -p msismall
#SBATCH -o /home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/clean_dscalars/output_logs/clean_dscalar_%j.out
#SBATCH -e /home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/clean_dscalars/output_logs/clean_dscalar_%j.err
#SBATCH -A btervocl

SUB=${1}
SES=${2}
TASK=${3}

# user specified inputs
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "job: ${SLURM_JOB_ID}"

dir_in=/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/derivatives/template_matching/
dscalar_infile=${dir_in}/sub-${SUB}/ses-${SES}/func/sub-${SUB}_ses-${SES}_task-${TASK}MENORDICtrimmed_space-fsLR_den-91k_desc-denoised_bold_spatially_interpolated_smoothed2mm_censoredFD2mm_zscoreRH_pythonTM.dscalar.nii
X="addpath(genpath('/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/clean_dscalars')); clean_dscalars_by_size('${dscalar_infile}',[],[],[],[],30,[],0,0,1);"
matlab_exec=/common/software/install/migrated/matlab/R2019a/bin/matlab
RandomHash=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16`
Tempmatlabcommand="matlab_command""$RandomHash"".m"

if [ -f "matlab_command""$RandomHash"".m" ]
then
	#echo "matlab_command.m found removing â¦"
	rm -fR "matlab_command""$RandomHash"".m"
fi

echo ${X} > "matlab_command""$RandomHash"".m"
cat "matlab_command""$RandomHash"".m"
${matlab_exec} -nodisplay -nosplash < "matlab_command""$RandomHash"".m"
rm -f "matlab_command""$RandomHash"".m"

echo "script finished at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "****************************************************"


