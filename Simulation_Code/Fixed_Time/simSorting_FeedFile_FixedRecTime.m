function [outputVars] = simSorting_FeedFile_FixedRecTime(inputFilePath,jobID,taskID,timeLimit)

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


%%

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

        % Calculate the variables
            % Find the cell type of each cell
                cellTypeTarget = cellType_onGrid(...
                    linIdx_ofPositions_onGrid_RepForNeighbors);
            % Find the cell types of the neighbors for each cell
                cellTypeNeighbor = cellType_onGrid(linIdx_ofNeighbors_onGrid);
            % Count the number of same cell type neighbors
                numSameCellType = accumarray(...
                    linIdx_ofPositions_onGrid_RepForNeighbors,...
                    cellTypeTarget==cellTypeNeighbor,[numCells,1]);

        % Pre-allocate space to save the output
            % Pre-allocate space for the cell IDs
                cellID_onGrid_OverTime = nan([sizeX,sizeY,numTP2Save]);
                cellID_onGrid_OverTime(:,:,1) = cellID_onGrid;
            % Pre-allocate space for the cell type
                cellType_onGrid_OverTime = nan([sizeX,sizeY,numTP2Save]);
                cellType_onGrid_OverTime(:,:,1) = cellType_onGrid;
            % Pre-allocate space for the number of same cell type
            % neighbors
                numSameCellTime_OverTime = nan([numCells,numTP2Save]);
                numSameCellTime_OverTime(:,1) = numSameCellType;
            % Pre-allocate space for the recorded time points
                timeInMin_OverTime = nan([numTP2Save,1]);
                timeInMin_OverTime(1,1) = 0;
            % Pre-allocate space for the time step size
                timeStepInMin_OverTime = nan([numTP2Save,1]);
            % Pre-allocate space for the total simulation time
                simTimeInSec_OverTime = nan([numTP2Save,1]);
                simTimeInSec_OverTime(1,1) = 0;


%% Run the simulation
    
    tic
    
    while timeSimInMin < time2SimInMin
    
        % Calculate the probability of swapping based on the number of neighbors that are the same cell type
        
            % Calculate the adhesion energy for each cell
                % Find the cell type of each cell
                    cellTypeTarget = cellType_onGrid(...
                        linIdx_ofPositions_onGrid_RepForNeighbors);
                % Find the cell types of the neighbors for each cell
                    cellTypeNeighbor = cellType_onGrid(linIdx_ofNeighbors_onGrid);
            
                % Count the number of same cell type neighbors
                    numSameCellType = accumarray(...
                        linIdx_ofPositions_onGrid_RepForNeighbors,...
                        cellTypeTarget==cellTypeNeighbor,[numCells,1]);
                    numDiffCellType = numNeighbors_ByPosition - numSameCellType;
        
                % Calculate the barrier height for each cell
                    E_Barrier = -((E_homo.*numSameCellType)+...
                        (E_het.*numDiffCellType));
    
            % Calculate the total energy barrier for the swap (sum of the adhesion
            % energy for both cells)
                E_Barrier_Swap = sum(E_Barrier(uniquePairs_linIdx_onGrid),2);
            % Calculate the rate of swapping for each pair
                rSwap = rSwap_0.*exp(-E_Barrier_Swap./kT_eff);
            % Set the initial timestep to capture an appropriate number of
            % events, set by pSwap_0
                dtInSec = min([pSwap_0/mean(rSwap),1/max(rSwap),time2SampleInMin*60,(nextTime2SampleInMin-timeSimInMin)*60]);
            % Calculate the probability of swapping
                pSwap = rSwap.*dtInSec;
    
        if any(pSwap>1)
            pSwap = pSwap/max(pSwap);
            dtInSec = dtInSec/max(pSwap);
        end
    
        if any(isnan(pSwap))
            ginput()
        end
    
        % For each pair of neighbors, calculate a random number between 0 and 1. 
        % If that number is less than the probability, they will switch. 
        % Because each pair is represented twice in the matrix, drop the switch
        % probability by a factor of two.
            randVals = rand([numPairs 1]);
            willSwitch = (randVals<pSwap);
    
        if any(willSwitch)
    
            % Find the linear indices for each cell in the pair        
                linIdx_Neighbor1_Swap = uniquePairs_linIdx_onGrid(willSwitch,1);
                linIdx_Neighbor2_Swap = uniquePairs_linIdx_onGrid(willSwitch,2);
    
            % Determine if any cell swaps more than once
                overlappingSwaps = (length([linIdx_Neighbor1_Swap ; linIdx_Neighbor2_Swap]) ~= ...
                    length(unique([linIdx_Neighbor1_Swap ; linIdx_Neighbor2_Swap]))); 
    
            % If there are cells chosen to swap more than once, reduce the
            % timestep by half and try again
            while overlappingSwaps
    
                % Drop the timestep by a factor of two
                    dtInSec = dtInSec/2;
    
                % Drop the probability by a factor of two
                    pSwap = pSwap./2;
    
                % For each pair of neighbors, calculate a random number between 0 and 1. 
                % If that number is less than the probability, they will switch. 
                % Because each pair is represented twice in the matrix, drop the switch
                % probability by a factor of two.
                    willSwitch = (randVals<pSwap);
                
                % Find the linear indices for each cell in the pair        
                    linIdx_Neighbor1_Swap = uniquePairs_linIdx_onGrid(willSwitch,1);
                    linIdx_Neighbor2_Swap = uniquePairs_linIdx_onGrid(willSwitch,2);
            
                % Determine if any cell still swaps more than once
                    overlappingSwaps = (length([linIdx_Neighbor1_Swap ; linIdx_Neighbor2_Swap]) ~= ...
                        length(unique([linIdx_Neighbor1_Swap ; linIdx_Neighbor2_Swap])));
        
            end
        
            % Find the cell IDs and cell types of the neighbor pair
                cellID_Neighbor1_Swap = cellID_onGrid(linIdx_Neighbor1_Swap);
                cellID_Neighbor2_Swap = cellID_onGrid(linIdx_Neighbor2_Swap);
                cellType_Neighbor1_Swap = cellType_onGrid(linIdx_Neighbor1_Swap);
                cellType_Neighbor2_Swap = cellType_onGrid(linIdx_Neighbor2_Swap);
        
            % Swap their identities
                cellID_onGrid(linIdx_Neighbor1_Swap) = cellID_Neighbor2_Swap;
                cellID_onGrid(linIdx_Neighbor2_Swap) = cellID_Neighbor1_Swap;
                cellType_onGrid(linIdx_Neighbor1_Swap) = cellType_Neighbor2_Swap;
                cellType_onGrid(linIdx_Neighbor2_Swap) = cellType_Neighbor1_Swap;
    
        end
    
        % Update the time and time step
            timeSimInMin = timeSimInMin + (dtInSec/60);
            tp = tp + 1;
    
        if timeSimInMin>=nextTime2SampleInMin
    
            % Record the time in minutes
                timeInMin_OverTime(tp2Record) = timeSimInMin;
            % Record the simulation time
                simTimeInSec_OverTime(tp2Record) = toc;
            % Record the time step in minutes for this timepoint
                timeStepInMin_OverTime(tp2Record) = dtInSec/60;
            % Record the sorting information
                cellID_onGrid_OverTime(:,:,tp2Record) = cellID_onGrid;
                cellType_onGrid_OverTime(:,:,tp2Record) = cellType_onGrid;
                numSameCellTime_OverTime(:,tp2Record) = numSameCellType;
            % Update the recording counter
                tp2Record = tp2Record + 1;
            % Update the next recording time in min
                if tp2Record<=length(times2SampleInMin)
                    nextTime2SampleInMin = times2SampleInMin(tp2Record);
                end
        end

        % Exit the simulation if the simulation time exceeds some value
        if toc >= timeLimitInSec
            ME3 = 'Time out';
            break
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
        end
    
end


