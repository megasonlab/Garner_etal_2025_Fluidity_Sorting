function [neighborsExchangedThisTP, slopes, times2Fit, window, slopeMedian, time_Min_NeighborExchange] = calculateNeighborExchangeRate(trackPositionsXYZ_UM_reshaped,time_Min_reshaped,neighborDistance_UM)


%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/19

% This function calculates the neighbor exchange rate in the simulations

%% Find the all nuclei within a fixed distance

disp('Finding neighbors')

% Calculate the number of track
    numTracks = size(trackPositionsXYZ_UM_reshaped,1);

% Pull out the time vector associated with the neighbor excange rate        
    time_Min_NeighborExchange = nansum(squeeze(time_Min_reshaped(:,1:(end-1))),1)./...
        sum(~isnan(squeeze(time_Min_reshaped(:,1:(end-1)))),1);
        
% Calculate the time step
    timeStep_min = diff(time_Min_reshaped,1,2);

% Calculate the number of timepoints
    maxNumTimepoints = size(trackPositionsXYZ_UM_reshaped,3);

% Preallocate space to store the number of neighbor exchanges
    neighborIdx_All = [];

% Find the nearest neighbobrs
    for numTimePoint = 1:maxNumTimepoints

        % Find the neighbors within the maximum distance
            neighborIdx = rangesearch(squeeze(trackPositionsXYZ_UM_reshaped(:,:,numTimePoint)),...
                squeeze(trackPositionsXYZ_UM_reshaped(:,:,numTimePoint)),...
                neighborDistance_UM(numTimePoint));

        % Add this to the cell array
            neighborIdx_All = [neighborIdx_All neighborIdx];

    end

% Count the number of neighbors
    numNeighbors = cellfun(@numel, neighborIdx_All);
    numNeighbors(numNeighbors==0) = nan;


%% Calculate the set differences

disp('Counting neighbor exchanges')

% Preallocate space to store the number of neighbor exchanges
    neighborExchangeRate = nan([numTracks,(maxNumTimepoints-1)]);
    neighborsExchangedThisTP = nan([numTracks,(maxNumTimepoints-1)]);

for numTrack = 1:numTracks

    % Pull out the timepoints for this track
        numTimePoints = find(squeeze(~isnan(time_Min_reshaped(numTrack,:))));
        numTimePoints = numTimePoints(2:end);

    for numTimePoint = numTimePoints

        % Record the number of neighbors exchanged
           neighborsExchangedThisTP(numTrack,numTimePoint-1) = ...
               length([setdiff(neighborIdx_All{numTrack,numTimePoint-1},...
               neighborIdx_All{numTrack,numTimePoint}) ...
               setdiff(neighborIdx_All{numTrack,numTimePoint},...
               neighborIdx_All{numTrack,numTimePoint-1})]);
    end


end

%% Perform moving linear regression


disp('Performing moving linear regression')

% Calculate the cumulative number of neighbors exchanged
    cumNumNeighborsExchanged = cumsum(neighborsExchangedThisTP,2);
% Calculate the cumulative number of neighbors exchanged
    meanCumNumNeighborsExchanged = mean(cumsum(neighborsExchangedThisTP,2),1,'omitnan');
    medianCumNumNeighborsExchanged = median(cumsum(neighborsExchangedThisTP,2),1,'omitnan');
% Choose the times to fit
    times2Fit = find((median(neighborsExchangedThisTP,1)<4)&(medianCumNumNeighborsExchanged>0));
% Choose the time windows to average over in # timepoints
    window = 15;

% Count the number of time points
    numTimePoints2Fit = length(times2Fit);

% Initialize vectors to store the slopes and interceptsot
    slopes = nan(numTracks,numTimePoints2Fit);
    intercepts = nan(numTracks,numTimePoints2Fit);

% Loop through each track and perform the moving linear regression
for cellNum = 1:numTracks
    [slopes(cellNum,:), intercepts(cellNum,:)] = ...
        movingLinearRegression(time_Min_NeighborExchange(times2Fit),...
        cumNumNeighborsExchanged(cellNum,times2Fit),window);
end

% Calculate the moving slope of the median
    [slopeMedian, interceptsMedian] = ...
        movingLinearRegression(time_Min_NeighborExchange(times2Fit),...
            meanCumNumNeighborsExchanged(times2Fit),window);
