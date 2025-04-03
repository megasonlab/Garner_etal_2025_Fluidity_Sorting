function [inputVars] = initializeFile(inputVars)      


%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by setUpAndRun_Sorting*.m, in order to initialize
% the simulation file
  
   % Create a template file to copy
        % Create a file to save the results and transfer to the cloud
        % Pull out the file name
            inputVars.fileNameOutTemp = [sprintf(inputVars.globalInfo.fileNameForSprintf,...
                inputVars.globalInfo.dateNum,inputVars.parameterValsNum,1) '_Temp'];
        % Concatenate to make the file path        
            inputVars.filePathOutTemp = [inputVars.globalInfo.folderPathOnLocal '/' inputVars.fileNameOutTemp];        
        % Save the file
            save(inputVars.filePathOutTemp,'-struct','inputVars','-v7.3')


end
