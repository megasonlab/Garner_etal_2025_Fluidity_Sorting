function [slopes, intercepts] = movingLinearRegression(xVals,yVals,window)


%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/19


% This function calculates a moving best fit slope and intercept in time
% series data


% Calculate the number of data points
    numPoints = length(xVals);

% Create a space to store the best fit values around each timepoint
    bestFitVals = cell(1, numPoints);

% Calculate the half window size
    halfWindowSize = floor(window/2);

% Loop through each point and calculate the local best fit values around each
% timepoint
for pointVal = 1:numPoints
    % Pull out the points to fit
        points2Fit = (pointVal-halfWindowSize):(pointVal+halfWindowSize);
    % Delete points outside of the vector
        points2Fit(points2Fit<=0) = [];
        points2Fit(points2Fit>numPoints) = [];
    % Perform linear regression
        bestFitVals{pointVal} = polyfit(xVals(points2Fit), yVals(points2Fit), 1)';
end

% Reformat the best fit values
    bestFitVals = cell2mat(bestFitVals)';
% Pull out just the slopes
    slopes = bestFitVals(:,1);
% Pull out the intercepts
    intercepts = bestFitVals(:,2);

end