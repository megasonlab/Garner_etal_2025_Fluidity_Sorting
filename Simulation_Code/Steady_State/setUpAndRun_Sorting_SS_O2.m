%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% The following script specifies the simulation parameters for the 
% parameter scan, initializes a simulation for each choice of parameters, 
% and then saves the initialized simulation variables to a
% file. The file path can then be input into the main simulation function
% simSorting_FeedFile_SS.m. You have the option to iterate 
% through a range of parameter values for one variable,
% as well as replicates for each choice of parameter value.

%%%%%% Run chooseModelParameters.m to set the base choice of parameters
% before running this code %%%%%


% Clear the system
    close all; 
    clear all;
    
% Choose a set of parameters values to scan through (you will assign 
% each parameter to its specific variable later on inside the for loop)

    % Choose the file path of default parameter vales
        baseParametersFilePath = ['\Garner_etal_2025_Fluidity_Sorting\'...
            'Simulation_Code\basicModelParameters.mat'];

    % Set the parameter values which will be updated in the for loop
        % Homotypic adhesion energy in multiples of laboratory kT
            parameterVals1 = 10.^[4, 5, 6, 7];
        % Viscosity multiplier (relative to wateR)
            parameterVals2 = [1, 10, 100, 1000];
        % Effective temperature in multiples of E_hommo
            parameterVals3 = fliplr([0.5, 0.75, 1, 1.5, 2]);
        % For now just choose one set
                parameterVals1 = 10.^[8];
                parameterVals2 = [1];
                parameterVals3 = linspace(1,2,21);
                        
    % Choose the number of replicates
        numReplicates = 5;

% Choose where to find/save variables

    % Choose the parent folder path to create the new folder in
        folderParentPathOnLocal = ['\Garner_etal_2025_Fluidity_Sorting\'...
            'Simulation_Data\Steady_State\'];
    
    % Create a file prefix to save the results (Change variable name "M" as
    % needed)
        fileNameForSprintf = 'Results_%s_SortSim_EhetEtaKT_%i_Rep_%i_in';            

% Set up the parameter scan based on the above user input

    % Create a list of all possible combinations of these parameter choices
        paramCombos = combvec(1:length(parameterVals1),...
            1:length(parameterVals2),...
            1:length(parameterVals3))';

    % Generature a list of the parameter combinations
        parameterVals = (1:size(paramCombos,1))';

    % Count the number of parameter sets
        numParameterSets = length(parameterVals);

    % Count the total number of simulations
        numSims = numParameterSets*numReplicates;

    % Create random seeds for each simulation
        % Shuffle the random number generator based on the current time
            rng('shuffle');
        % Randomly assign seeds
            seedVals = nan(numParameterSets,numReplicates);
            seedVals(:) = randsample(numSims,numSims);
        
    % Create a folder to save the results
        % Pull out the date and convert to string for folder naming convention
            str = date();
            dateNum = datestr(str,'yyyy/mm/dd');
            dateNum = dateNum(dateNum~='/');
        % Create the new folder name and path, assuming it doesn't already
        % exist, otherwise append a number to the end of the date until we have
        % a new folder
            folderNameOnLocal = sprintf('Results%s',dateNum);
            folderPathOnLocal = [folderParentPathOnLocal folderNameOnLocal];
            if ~exist(folderPathOnLocal,'dir')
                mkdir(folderPathOnLocal);
            else
                n=1;
                while exist(folderPathOnLocal,'dir')
                    folderNameOnLocal = sprintf('Results%s_%i',dateNum,n);
                    folderPathOnLocal = [folderParentPathOnLocal folderNameOnLocal];
                    n=n+1;
                end
                mkdir(folderPathOnLocal);
            end
    
        % Record this file save location (the strange nomenclature is for more 
        % general applications wusing parallelization)
            folderPathOnCloud = folderPathOnLocal;
            destinationFolderPath = folderPathOnCloud;
        
    % Put this info into a structure
        % Make the structure
            globalInfo = v2struct;
        % Clear the remaining variables
            clearvars -except globalInfo
        
% Submit the jobs        
for parameterValsNum = 1:length(globalInfo.parameterVals)

    if parameterValsNum==1
        % Copy the functions necessary to run this script into the folder
            copyfile(which('v2struct.m'),globalInfo.folderPathOnLocal)
            copyfile(which('simSorting_FeedFile_SS.m'),globalInfo.folderPathOnLocal)     
            copyfile(which('createFitAsymHill.m'),globalInfo.folderPathOnLocal)      
            copyfile(which('evalInverseAsymHill.m'),globalInfo.folderPathOnLocal)     
    end
    
    % Set up the simulation for this parameter set (excluding any variables
    % that need to be chosen randomly, like the initial filament positions)
    
        % Load the base set of parameters
            load(globalInfo.baseParametersFilePath)

        % Adjust whichever parameters we're scanning through
            % Homotypic energy (in # kT_lab)
                E_homo = -globalInfo.parameterVals1(globalInfo.paramCombos(parameterValsNum,1));
            % Effective kT as a fraction of E_homo (in # kT_lab)
                kT_eff = -globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3))*E_homo;
            % Choose the viscosity of the tissue in pN nm ms
                etaVal = globalInfo.parameterVals2(globalInfo.paramCombos(parameterValsNum,2))*10^(-3);

        % Load these variables into a structure
            inputVars = v2struct();
            % Clear all other variables
                clearvars -except globalInfo parameterValsNum inputVars

        % Initialize the simulation (must be in this order)
            % Initialize the biophysical parameters
                inputVars = initializeBiophysicalParameters(inputVars);
            % Initialize the cells on the grid
                inputVars = initializeGrid(inputVars);
            % Initialize the recording information
                inputVars = initializeRecording(inputVars);
            % Pre-allocate space for large variables in a file    
                inputVars = initializeFile(inputVars);

    % Submit this function as a batch job for each replicate
        for replicateNum = 1:globalInfo.numReplicates

            % Print which loop we're on
                sprintf('Submitting replicate %i of %i for parameter set %i of %i',...
                    replicateNum,globalInfo.numReplicates,parameterValsNum,length(globalInfo.parameterVals))

            % Copy the file for this replicate
                inputVars = copyFileForReplicate(inputVars,globalInfo,parameterValsNum,replicateNum);

        end

        delete([inputVars.filePathOutTemp '.mat'])

    % Clear the remaining variables
        clearvars -except globalInfo parameterValsNum

end
