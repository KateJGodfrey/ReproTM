#!/bin/bash -l

#SBATCH -J matlab
#SBATCH --ntasks=4
#SBATCH --tmp=40gb
#SBATCH --mem=40gb
#SBATCH -t 1:00:00
#SBATCH -p msismall
#SBATCH -o output_logs/matlabTM_%j.out
#SBATCH -e output_logs/matlabTM_%j.err
#SBATCH -A btervocl

# overview 
# this shell script calls a hardcoded version of the original matlab template matching
# this script is potentially useful for comparision between matlab and python implimentations

# user specified inputs
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "job: ${SLURM_JOB_ID}"


X="addpath(genpath('/home/btervocl/shared/projects/MINT/kgodfrey/template_matching')); template_matching_hardcode_RH"
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


