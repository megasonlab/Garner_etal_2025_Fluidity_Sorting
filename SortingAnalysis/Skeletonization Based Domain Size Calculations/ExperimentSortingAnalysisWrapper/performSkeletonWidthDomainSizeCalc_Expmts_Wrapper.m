
%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This script takes a folder containing experimental segmented image data,
% finds each of the segmented image files, and for each file calls
% performSkeletonWidthDomainSizeCalc_Expmts_FeedFile_O2 to perform the
% domain size calculation

% Clear the system
    close all;
    clear all;

% Input the image file path
    dataFolderPaths = ['\Garner_etal_2025_Fluidity_Sorting\Experimental_Data\Processed_Data\240418_spheroid_assay\'];
   % dataFolderPaths = ['\Garner_etal_2025_Fluidity_Sorting\Experimental_Data\Processed_Data\240420_spheroid_assay\'];

% Initialize the counters
    numCompleted = 0;
    numFailed = 0;
    numNeverFinished = 0;

for folderNum = 1:length(dataFolderPaths)

    % Print the folder number
        folderNum

    % Pull out the folder
        dataFolderPath = dataFolderPaths{folderNum};
    % Pull out the mat files
        matFiles = dir([dataFolderPath '*_MaxIP_Segmentation.mat']);
    % Sort the files the way a human would
       [matFiles] = natsortfiles(matFiles);

    for matFileNum = 1:length(matFiles)

        % Print the mat file num
            matFileNum
    
        % Pull out the filename
           inputFilePath = [dataFolderPath matFiles(matFileNum).name];
    
        % Load the result
            clear domainWidthinMicrons_RepByPixel_CellType2
            load(inputFilePath,'domainWidthinMicrons_RepByPixel_CellType2')

        if exist('domainWidthinMicrons_RepByPixel_CellType2','var') == 1

            % Update the counter
                numCompleted = numCompleted+1;

        else

            clear ME20
            load(inputFilePath,'ME20')

            if exist('ME20','var') == 1

                % Update the counter
                    numFailed = numFailed+1;

                disp('Analysis failed')

            else

                % Update the counter
                    numNeverFinished = numNeverFinished+1;

                disp('Analysis never completed')
        
                % Run the analysis
                tic
                    performSkeletonWidthDomainSizeCalc_Expmts_FeedFile_O2(inputFilePath);
                toc
    
            end
        end

   end

end

% Initialize the counters
    numCompleted
    numFailed
    numNeverFinished


