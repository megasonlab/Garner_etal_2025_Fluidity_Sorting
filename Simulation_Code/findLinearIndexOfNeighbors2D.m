function [linIdx_ofPositions_onGrid_RepForNeighbors, linIdx_ofNeighbors_onGrid] = findLinearIndexOfNeighbors2D(sizeGrid,connectivity,BCs)


%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by initializeGrid.m, in order to initialize
% the linear indices of each cell's neighbors


% sizeGrid
% Connectivity: Connectivity can be 4-connected (only cells who share faces)
% or 8-connected (cells sharing either a face or a corner)
% BCs: Boundary condition can be periodic (0) or fixed (1) 


% Find the linear index of the neighbors for each pixel on the grid
    % Pull out the size of the grd
        sizeX = sizeGrid(1);
        sizeY = sizeGrid(2);
    % Calculate the number of cells
        numCells = sizeX*sizeY;
    % Create a linear index for each position on the grid
        linIdx_ofPositions_onGrid = 1:numCells;
    % Identify the subscript indices for each matrix entry
        [row,col] = ind2sub([sizeX, sizeY],linIdx_ofPositions_onGrid);
    % Identify the [row col] subscript indices of the neighbors 
    % (relative to the pixel of interest)
        if connectivity==4
            subIdx_Neighbor_Rel = [-1 0; 0 -1; 0 1; 1 0];
        elseif connectivity==8
            subIdx_Neighbor_Rel = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 0 1];
        else
            'Connectivity must be 4-connected (4) or 8-connected (8) '
            return
        end
    % Calculate the number of neighbors
        maxNumNeighbors = size(subIdx_Neighbor_Rel,1);
    % Create a repeated matrix of the cell IDs for each neighbor
        linIdx_ofPositions_onGrid_RepForNeighbors = repelem(linIdx_ofPositions_onGrid(:),maxNumNeighbors);
    % Find the subscript indices of the neighbors
        % Make a repeated relative neighbor index matrix for each pixel in
        % the matrix
            subIdx_Neighbors_Rel_RepAllPixels = repmat(subIdx_Neighbor_Rel,[numCells,1]);
        % Make a repeated index matrix for each neighbor for each pixel in
        % the matrix
            subIdx_Pixels_RepForNeighbors = [repelem(row,maxNumNeighbors)', repelem(col,maxNumNeighbors)'];
        % Add the subscript indices of the pixel position with the relative
        % neihgbor position together to find the absolute neighbor subscript 
        % indices
            subIdx_Neighbors_ForPixels = subIdx_Pixels_RepForNeighbors + subIdx_Neighbors_Rel_RepAllPixels;
        % Implement boundary conditions
        if BCs==0
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,1)<1,1) = ...
                subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,1)<1,1) + sizeX;
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,1)>sizeX,1) = ...
                subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,1)>sizeX,1) - sizeX;
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,2)<1,2) = ...
                subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,2)<1,2) + sizeY;
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,2)>sizeY,2) = ...
                subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,2)>sizeY,2) - sizeY; 
        elseif BCs==1            
            linIdx_ofPositions_onGrid_RepForNeighbors(subIdx_Neighbors_ForPixels(:,1)<1) = [];
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,1)<1,:) = [];
            linIdx_ofPositions_onGrid_RepForNeighbors(subIdx_Neighbors_ForPixels(:,1)>sizeX) = [];
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,1)>sizeX,:) = [];
            linIdx_ofPositions_onGrid_RepForNeighbors(subIdx_Neighbors_ForPixels(:,2)<1) = [];
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,2)<1,:) = [];      
            linIdx_ofPositions_onGrid_RepForNeighbors(subIdx_Neighbors_ForPixels(:,2)>sizeY) = []; 
            subIdx_Neighbors_ForPixels(subIdx_Neighbors_ForPixels(:,2)>sizeY,:) = []; 
        else
            'Boundary conditions must be periodic (0) or fixed (1) '
            return
        end  
    % Convert back to linear indices
        linIdx_ofNeighbors_onGrid = sub2ind([sizeX, sizeY],subIdx_Neighbors_ForPixels(:,1),subIdx_Neighbors_ForPixels(:,2));
    % Calculate the number of pairs
        numPairs = length(linIdx_ofNeighbors_onGrid);

end