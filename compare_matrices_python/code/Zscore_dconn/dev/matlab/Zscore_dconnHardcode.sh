#!/bin/bash -l

#SBATCH -J matlabZ
#SBATCH --ntasks=4
#SBATCH --tmp=120gb
#SBATCH --mem=120gb
#SBATCH -t 1:00:00
#SBATCH -p msibigmem
#SBATCH -o output_logs/matlabZ_%j.out
#SBATCH -e output_logs/matlabZ_%j.err
#SBATCH -A btervocl

# user specified inputs
echo "****************************************************"
echo "script started at: $(date +"%Y-%m-%d") $(date +"%r")"
echo "job: ${SLURM_JOB_ID}"


X="addpath(genpath('/home/btervocl/shared/projects/MINT/kgodfrey/template_matching/code/zscore_dconn')); Zscore_dconnHardcode"
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


