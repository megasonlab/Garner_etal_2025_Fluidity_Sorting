function [radial_Average, radial_STD] = calculateRadialAverage(matrix,centerPositionInPixels)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This function performs a radial average of a matrix, based on a user-specified center position
% Inputs: matrix - a 2D matrix 
%         centerPositionInPixels - the center position in pixels 
%             following the format [rowNum colNum]
% Outputs: radial_Average - a 1 x max(rowNum,colNum) vector containing the
%           radial average as a function of # pixels away from the center


if ~exist('centerPositionInPixels','var')
	centerPositionInPixels = [1 1];
end

% For each pixel in the matrix, calculate the distance from the 
% specified center (in number of pixels)
    % Calculate the matrix size in pixels [numRows numCols]
        matrixSizeInPixels = size(matrix);
    % Create the matrix of row and column positions in pixels
        [Col_Pos_Pix,Row_Pos_Pix] = meshgrid(1:matrixSizeInPixels(2),1:matrixSizeInPixels(1));    
    % Calculate the radius of each pixel from the center position
        R_Pix = sqrt((Row_Pos_Pix-centerPositionInPixels(1)).^2 + ...
            (Col_Pos_Pix-centerPositionInPixels(2)).^2);

% For each pixel in the matrix, calculate the distance from the 
% specified center in discretized radial steps
    % Choose the step size for the radii to check in pixels
        radial_Step = 1;
    % For each pixel in the matrix, calculate the distance from the 
    % centroid (in number of radial steps)
        radius_Discrete = round(R_Pix/radial_Step)+1;
    % Calculate the maximum radial distance in pixels
        max_Radius_Pixels = max(radius_Discrete(:));

% Calculate the radial average of the matrix for these discretized radial steps
    radial_Average = accumarray(radius_Discrete(:),matrix(:),...
        [max_Radius_Pixels 1],@mean,nan);
    radial_STD = accumarray(radius_Discrete(:),matrix(:),...
        [max_Radius_Pixels 1],@std,nan);

