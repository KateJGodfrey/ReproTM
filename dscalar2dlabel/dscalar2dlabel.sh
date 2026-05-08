#!/bin/bash

# simple script to convert dscalar to dlabel with workbench

# input dscalar
in_dscalar=${1}
# output dlabel
out_dlabel=${2}
# path to label file
label_file=*/tpl-ABCC2026-a3-9to16_15networks_label_list.txt
# path to connectome workbench wb_command
wb_command=/path/connectome/workbench/wb_command

# run wb_command to convert dscalar to dlabel
${wb_command} -cifti-label-import ${in_dscalar} ${label_file} ${out_dlabel} 