#!/bin/bash

# Declare an array variable named filePaths
declare -a filePaths		

# Fill in the array from the input arguments
filePaths=( "$@" )			

# Load the matlab module for a specified version
module load matlab/2023a		

# Extract the folder path to the data
folderPath="$(dirname "${filePaths[$((SLURM_ARRAY_TASK_ID - 1))]}")"
					 
# Load matlab, set the matlab path to the current folder, and run the function using the file path, job ID, and task ID as input arguments
matlab -batch "addpath(\"$folderPath\",'-begin'); performSkeletonWidthDomainSizeCalc_Sims_FeedFile_O2(\"${filePaths[$((SLURM_ARRAY_TASK_ID - 1))]}\")"

					
