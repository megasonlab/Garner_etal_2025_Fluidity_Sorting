function [outputVars] = simSorting_FeedFile_Ising_LocalSwap(inputFilePath,jobID,taskID,timeLimit)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% The following code performs a 2D grid-based simulation of adhesion-based sorting

% Inputs: inputFilePath - the file path containing the initialized simulation 
%           jobID - the O2 SLURM Job ID 
%           taskID - the O2 SLURM Task ID 
%           timeLimit - time limit of slurm job in D-HH_MM_SS format 

% Outputs: outputVars - a structure containing the input parameters as well
% as the simulation results. outputVars is saved directly to the hard drive, 
% using the save file path specified in inputVars, and then the structure 
% is deleted (in the MATLAB workspace) if the save is performed without
% errors.

% This code was written by Rikki M. Garner and was last updated 2024/10/13

%%
tic
% Create an empty error message
    ME=[];
    ME3=[];

% Try to run the function, otherwise print the error message
try

%% Set up the simulation

    % Convert the filepath from string to char
        inputFilePath = convertStringsToChars(inputFilePath);

    % Choose the output file path
        outputFilePath = sprintf('%s_%08d_%08d_out.mat',inputFilePath(1:(end-7)),jobID,taskID);

    % Parse the time limit
        % Split up the string using delimiters
            timeLimitSplit = strsplit(timeLimit,{'-',':'});
        % Convert to numerical values
            timeLimitSplit = str2double(timeLimitSplit);
        % Convert to seconds
            if length(timeLimitSplit)==2
                timeLimitInSec = (timeLimitSplit(1)*60) + timeLimitSplit(2);
            elseif length(timeLimitSplit)==3
                timeLimitInSec = (timeLimitSplit(1)*60*60) + (timeLimitSplit(2)*60) + timeLimitSplit(3);
            elseif length(timeLimitSplit)==4
                timeLimitInSec = (timeLimitSplit(1)*24*60*60) + (timeLimitSplit(2)*60*60) + (timeLimitSplit(3)*60) + timeLimitSplit(4);
            end
        % Leave 5 minutes to save the data
            timeLimitInSec = max(timeLimitInSec*0.95,timeLimitInSec-(5*60));

    % Load the variables from the specified file
        load(inputFilePath)

    % Reset the random number generator seed
        rng(globalInfo.seedVals(parameterValsNum,replicateNum))
        seedInfo = rng;

    % Initialize the simulation

        % Assign exactly half of the cells, randomly sampled, to be one
        % cell type. The remaining cells will be the 2nd cell type.
            cellType_onGrid(:)= nan;
            cellType_onGrid(randsample(numCells,round(numCells/2))) = 1;
            cellType_onGrid(isnan(cellType_onGrid)) = 2;

        % Pre-allocate space to save the output
            % Pre-allocate space for the cell IDs
                cellID_onGrid_OverTime = nan([sizeX,sizeY,numTP2Save]);
                cellID_onGrid_OverTime(:,:,1) = cellID_onGrid;
            % Pre-allocate space for the cell type
                cellType_onGrid_OverTime = nan([sizeX,sizeY,numTP2Save]);
                cellType_onGrid_OverTime(:,:,1) = cellType_onGrid;
            % Pre-allocate space for the recorded time points
                timeInSteps_OverTime = nan([numTP2Save,1]);
                timeInSteps_OverTime(1,1) = 0;
            % Pre-allocate space for the total simulation time
                simTimeInSec_OverTime = nan([numTP2Save,1]);
                simTimeInSec_OverTime(1,1) = 0;

        % Pre-allocate the row and column indices
            [target_rows, target_cols]  = ind2sub(size(cellType_onGrid), 1:numel(cellType_onGrid));
        % Initialize the neighbors
            subIdx_NNShift = [0 -1; 1 0; 0 1; -1 0];

        % Calculate the number of time steps
            numTimeSteps = 100000*numel(cellType_onGrid);
            tp2Record = 2;

    % Choose the timepoints to sample
        % % Linear
           % times2SampleInMin = 0:(time2SimInMin/(numTP2Save-1)):time2SimInMin;
        % Log-distributed
            firstTime2Sample = 1;
            lastTime2Sample = numTimeSteps;
            interval = (log10(lastTime2Sample)-log10(firstTime2Sample))/(numTP2Save-2);
            timeSteps2Sample = unique([0 round(10.^(log10(firstTime2Sample):interval:log10(lastTime2Sample)))]);

%% Run the simulation
    
tic
    
for stepNum = 1:numTimeSteps

    %% Pick a random cell
        target_linIdx = randi(numel(cellType_onGrid));
    % Find its row and column indices
        target_row = target_rows(target_linIdx);
        target_col = target_cols(target_linIdx);
    % Pick a random neighbor
        NN_row_col = subIdx_NNShift(randi(4),:);
    % Find its row and column indices (assuming periodic BCs)
        NN_row = mod(target_row + NN_row_col(1) - 1, size(cellType_onGrid,1)) + 1;
        NN_col = mod(target_col + NN_row_col(2) - 1, size(cellType_onGrid,2)) + 1;

    % Calculate the change in energy due to the switch
        if cellType_onGrid(target_row, target_col)==cellType_onGrid(NN_row, NN_col)
            % If the cell types are the same, the change in energy is 0
            dE = 0;
        else
            % Find all nearest neighbors of the target cell
                target_above = mod(target_row - 1 - 1, size(cellType_onGrid,1)) + 1;
                target_below = mod(target_row + 1 - 1, size(cellType_onGrid,1)) + 1;
                target_left  = mod(target_col - 1 - 1, size(cellType_onGrid,2)) + 1;
                target_right = mod(target_col + 1 - 1, size(cellType_onGrid,2)) + 1;
            % Find all nearest neighbors of the neighbor cell
                NN_above = mod(NN_row - 1 - 1, size(cellType_onGrid,1)) + 1;
                NN_below = mod(NN_row + 1 - 1, size(cellType_onGrid,1)) + 1;
                NN_left  = mod(NN_col - 1 - 1, size(cellType_onGrid,2)) + 1;
                NN_right = mod(NN_col + 1 - 1, size(cellType_onGrid,2)) + 1;
            % Determine the cell types of the target cell, NN cell, and their neighbors
                target_cellType = cellType_onGrid(target_row, target_col);
                target_neighborCellTypes = [cellType_onGrid(target_above,target_col); ...
                    cellType_onGrid(target_row,target_left); ...
                    cellType_onGrid(target_row,target_right); ...
                    cellType_onGrid(target_below,target_col)];
                NN_cellType = cellType_onGrid(NN_row, NN_col);
                NN_neighborCellTypes = [cellType_onGrid(NN_above,NN_col); ...
                    cellType_onGrid(NN_row,NN_left); ...
                    cellType_onGrid(NN_row,NN_right); ...
                    cellType_onGrid(NN_below,NN_col)];
            % Calculate the number of neighbors of the same cell type
                numSameCellType_Current = sum(target_cellType==target_neighborCellTypes) + ...
                    sum(NN_cellType==NN_neighborCellTypes);
            % Calculate the change in energy if they were to switch
                dE = (6-(2*numSameCellType_Current))*(E_homo-E_het);
        end

        % Determine whether the cells will swap based on the difference in
        % energy
            % Calculate the probability of swapping
                pSwap = exp(-dE./kT_eff);
            if rand<pSwap
                % Find the cell IDs and cell types of the neighbor pair
                    cellID_Neighbor1_Swap = cellID_onGrid(target_row,target_col);
                    cellID_Neighbor2_Swap = cellID_onGrid(NN_row, NN_col);
                    cellType_Neighbor1_Swap = cellType_onGrid(target_row,target_col);
                    cellType_Neighbor2_Swap = cellType_onGrid(NN_row, NN_col);
                % Swap their identities
                    cellID_onGrid(target_row,target_col) = cellID_Neighbor2_Swap;
                    cellID_onGrid(NN_row, NN_col) = cellID_Neighbor1_Swap;
                    cellType_onGrid(target_row,target_col) = cellType_Neighbor2_Swap;
                    cellType_onGrid(NN_row, NN_col) = cellType_Neighbor1_Swap;
            end
   
    % Record intermittenly
   if stepNum==timeSteps2Sample(tp2Record)
   
        % Record the time in minutes
            timeInSteps_OverTime(tp2Record) = stepNum;
        % Record the simulation time
            simTimeInSec_OverTime(tp2Record) = toc;
        % Record the sorting information
            cellID_onGrid_OverTime(:,:,tp2Record) = cellID_onGrid;
            cellType_onGrid_OverTime(:,:,tp2Record) = cellType_onGrid;
        % Update the recording counter
            tp2Record = tp2Record + 1;

   end

end

catch errorMessage
    ME = errorMessage;
end


%% Save the results

    % Pack the data into a structure
        outputVars = v2struct;
    % Clear the remaining variables
        clearvars -except outputVars
    % Save the results
        try
            if isempty(outputVars.ME)&isempty(outputVars.ME3)
                % Save the data to a new output file
                    save(outputVars.outputFilePath,'-struct','outputVars','-v7.3');
                % Clear outputVars
                   % clear outputVars
            elseif isempty(outputVars.ME3)
                % Choose the output file path
                    outputVars.outputFilePath = [outputVars.outputFilePath(1:(end-4)) '_OtherFail.mat'];
                % Save the data to a new output file
                    save(outputVars.outputFilePath,'-struct','outputVars','-v7.3');
                % Clear outputVars
                  %  clear outputVars
            else
                % Choose the output file path
                    outputVars.outputFilePath = [outputVars.outputFilePath(1:(end-4)) '_TimeOutFail.mat'];
                % Save the data to a new output file
                    save(outputVars.outputFilePath,'-struct','outputVars','-v7.3');
                % Clear outputVars
                   % clear outputVars
            end
        catch ME2
            % Save the error
                outputVars.ME2 = ME2;
                ME2
        end
toc
end


