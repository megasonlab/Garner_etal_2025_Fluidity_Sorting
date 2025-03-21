function [ME20] = performSkeletonWidthDomainSizeCalc_Sims_FeedFile_O2(inputFilePath)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This function loads a simulated cell type matrix, calls 
% calculateDomainSize_SkelDist_SizeDeptFill to calulcate the 
% domain size mean and variance over time, and saves the results to the disk

% Create an empty error message
    ME20=[];

% Try to run the function, otherwise print the error message
try

    % Convert the filepath from string to char
        inputFilePath = convertStringsToChars(inputFilePath);

    % Load any existing analysis files
        clear domainWidthinMicrons_RepByPixel
        load(inputFilePath,'domainWidthinMicrons_RepByPixel')

if exist('domainWidthinMicrons_RepByPixel','var') == 1

else

    % Try to load any existing error messages
    clear ME20
    load(inputFilePath,'ME20')

    if exist('ME20','var') == 1

    else

    % Load the relevant variables
        load(inputFilePath,'cellType_onGrid_OverTime')
        load(inputFilePath,'timeInMin_OverTime')
        load(inputFilePath,'cellRadius')

    % Create a structure to store the data
        domainWidthinMicrons_RepByPixel = struct;
    % Choose the timepoints to analyze
        % All timepoints
            domainWidthinMicrons_RepByPixel.tp2Analze = 1:size(cellType_onGrid_OverTime,3);
        % Every 100th timepoint
            %domainWidthinMicrons_RepByPixel.tp2Analze = 1:100:size(cellType_onGrid_OverTime,3);
        % Every 30 minutes
            % % Calculate the target times
            %     timeInMin_Target = 0:30:max(timeInMin_OverTime);
            % % Calculate the time different between all possible pairs of
            % % the target time and the actual timepoints
            %     timeDifferenceAllPairs = abs(timeInMin_OverTime(:)-timeInMin_Target(:).');
            % % Find the closest simulated timepoint to each experimental
            % % timepoint
            %     [M, I] = min(timeDifferenceAllPairs,[],1);
            % % Keep only the unique timepoints
            %     domainWidthinMicrons_RepByPixel.tp2Analze = unique(I');
    % Choose the cell diameter
        domainWidthinMicrons_RepByPixel.cellDiameterInUM = cellRadius*2/1000;    
    % Input the target spatial resolution
        domainWidthinMicrons_RepByPixel.micronsPerPixel_X_Target = 1.2430;
    % Input the current spatial resolution
        domainWidthinMicrons_RepByPixel.micronsPerPixel_X_Original = domainWidthinMicrons_RepByPixel.cellDiameterInUM;
    % Choose whether to plot the data
        doPlot = false;

    % Run the fitting algorithm
        [domainWidthinMicrons_RepByPixel.Mean, domainWidthinMicrons_RepByPixel.STD, domainWidthinMicrons_RepByPixel.SE, ...
            domainWidthinMicrons_RepByPixel.Skew, domainWidthinMicrons_RepByPixel.Kurtosis] = ...
                calculateDomainSize_SkelDist_SizeDeptFill(cellType_onGrid_OverTime(:,:,domainWidthinMicrons_RepByPixel.tp2Analze), ...
                domainWidthinMicrons_RepByPixel.cellDiameterInUM, domainWidthinMicrons_RepByPixel.micronsPerPixel_X_Original, ...
                domainWidthinMicrons_RepByPixel.micronsPerPixel_X_Target, doPlot);

    % Save the analyzed data to the output file
        save(inputFilePath,'domainWidthinMicrons_RepByPixel','-append');

    end

end

catch errorMessage
    % Record the error message
        ME20 = errorMessage;
    % Save the error message to the output file
        save(inputFilePath,'ME20','-append');
end

end
