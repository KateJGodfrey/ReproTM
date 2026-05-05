#! /bin/sh

# enter code directory
#cd "$( dirname "${BASH_SOURCE[0]}" )"

## Matlab command and usage

# Zscore_dconn(input_cifti_name,output_cifti_name,path_wb_c)
X="addpath('/home/faird/shared/code/internal/utilities/Zscore_dconn'); Zscore_dconn('${1}', '${2}', '${3}')"

#Hermosillo R. 11/02/2023

###########################################################################################
matlab_exec=/common/software/install/migrated/matlab/R2019a/bin/matlab
RandomHash=`cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16`
Tempmatlabcommand="matlab_command""$RandomHash"".m"

if [ -f "matlab_command""$RandomHash"".m" ]
then
	#echo "matlab_command.m found removing â¦"
	rm -fR "matlab_command""$RandomHash"".m"
fi

#echo ${X} 
echo ${X} > "matlab_command""$RandomHash"".m"
cat "matlab_command""$RandomHash"".m"
${matlab_exec} -nodisplay -nosplash < "matlab_command""$RandomHash"".m"
rm -f "matlab_command""$RandomHash"".m"
