
%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This script finds all MATLAB files included in a directory and returns
% the list of ToolBoxes required to run the code


% Choose the directory to search for MATLAB files
    rootdir = '\Garner_etal_2025_Fluidity_Sorting\Code\NeighborAnalysis\';
% Make a list of all the MATLAB files in that directory and subdirectories
    filelist = dir(fullfile(rootdir, '**\*.m')); 
% Determine which toolboxes are required to run these MATLAB files
    [fList,pList] = matlab.codetools.requiredFilesAndProducts({filelist.name});
    {pList.Name}


% Choose the directory to search for MATLAB files    
    rootdir = '\Garner_etal_2025_Fluidity_Sorting\Code\SortingAnalysis\';
% Make a list of all the MATLAB files in that directory and subdirectories
    filelist = dir(fullfile(rootdir, '**\*.m')); 
% Determine which toolboxes are required to run these MATLAB files
    [fList,pList] = matlab.codetools.requiredFilesAndProducts({filelist.name});
    {pList.Name}


% Choose the directory to search for MATLAB files
    rootdir = '\Garner_etal_2025_Fluidity_Sorting\Code\Figure_Panels_Code\';
% Make a list of all the MATLAB files in that directory and subdirectories
    filelist = dir(fullfile(rootdir, '**\*.m')); 
% Determine which toolboxes are required to run these MATLAB files
    [fList,pList] = matlab.codetools.requiredFilesAndProducts({filelist.name});
    {pList.Name}

