#!/bin/bash

# Choose the folder contianing the files to run
folderPath="/home/rig121/Sorting_Fluidity/Parameter_Scans/Fixed_Time/Results20240820" 

# Choose the pattern for the files to run
pattern="*_out.mat"					

# Pull out all the files in the specified folder matching the specified pattern
filePaths=("$folderPath"/$pattern)			

# Count the number of files (to set the sbatch array limits)
numFiles="${#filePaths[@]}"			

# Loop through each file and submit the run as a batch job using the slurm scheduler
sbatch --array=1-"$numFiles" -c 1 -t 0-02:00 -p short --mem=1200M -o hostname_%A_%a.out -e hostname_%A_%a.err "${folderPath}/run_CalculateDomainSize_SkelDist_20240820.sh" "${filePaths[@]}"
#sbatch --array=1-1 -c 1 -t 0-00:35 -p short --mem=2000M -o hostname_%A_%a.out -e hostname_%A_%a.err "${folderPath}/run_CalculateDomainSize_SkelDist_20240221.sh" "${filePaths[@]}"

