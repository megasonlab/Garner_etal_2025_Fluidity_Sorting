function [inputVars] = copyFileForReplicate(inputVars,globalInfo,parameterValsNum,replicateNum)


%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by setUpAndRun_Sorting*.m, in order to copy
% replicate files

    % Create a file to save the results and transfer to the cloud
        % Pull out the file name
            fileNameOut = sprintf(globalInfo.fileNameForSprintf,...
                globalInfo.dateNum,parameterValsNum,replicateNum);
        % Concatenate to make the file path        
            filePathOutOnLocal = [globalInfo.folderPathOnLocal '/' fileNameOut];
            filePathOutOnCloud = [globalInfo.folderPathOnCloud '/' fileNameOut];
        % Add the file names to the inputParams structure
            inputVars.replicateNum = replicateNum;
            inputVars.fileNameOut = fileNameOut;
            inputVars.filePathOutOnLocal = filePathOutOnLocal;
            inputVars.filePathOutOnCloud = filePathOutOnCloud;
        % Save the file
            copyfile([inputVars.filePathOutTemp '.mat'], [filePathOutOnLocal '.mat'])
            save(filePathOutOnLocal,'-append','-struct','inputVars')

end