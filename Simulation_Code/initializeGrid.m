function [inputVars] = initializeGrid(inputVars)


%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by setUpAndRun_Sorting*.m, in order to initialize
% the cell grid

    % Unpack this structure and clear the original structure
        v2struct(inputVars)
        clear inputVars        
            
    % Initalize the grid
        % Calculate the number of cells
            numCells = sizeX*sizeY;
        % Assign each cell on the grid an ID
            cellID_onGrid = nan([sizeX,sizeY]);
            cellID_onGrid(:) = 1:numCells;
        % Assign each cell a cell type
            cellType_onGrid = nan([sizeX,sizeY]);

    % Create the linear indices for target and neighbor cell comparisons
        % Input the size of the grid
            sizeGrid = [sizeX sizeY];        
        % Calculate the linear index pairs for the targets and their neighbors
            [linIdx_ofPositions_onGrid_RepForNeighbors, linIdx_ofNeighbors_onGrid] = ...
                findLinearIndexOfNeighbors2D(sizeGrid,connectivity,BCs);
        % Calculate the number of pairs
            numTotalPairs = length(linIdx_ofPositions_onGrid_RepForNeighbors);
        % Calculate the number of neighbors for each 
            numNeighbors_ByPosition = ...
                accumarray(linIdx_ofPositions_onGrid_RepForNeighbors,1);
        % Find all possible pairs of neighbors (without replacement)
            % Sort the neighbors pairs such that the lower index neighbor is first
                sorted_LinIdx = sort([linIdx_ofPositions_onGrid_RepForNeighbors, linIdx_ofNeighbors_onGrid],2);
            % Find linear indices for the unique pairs
                [uniquePairs_linIdx_onGrid,I,J] = unique(sorted_LinIdx, 'rows', 'first');  
            % Count the number of pairs
                numPairs = size(uniquePairs_linIdx_onGrid,1);
    
    % Load these variables into a structure
        inputVars = v2struct();

end
