function [ME20] = performSkeletonWidthDomainSizeCalc_Expmts_FeedFile_O2(inputFilePath)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This function loads an experimental segmented image file, and alls 
% calculateDomainSize_SkelDist_SizeDeptFill to calulcate the 
% domain size mean and variance over time, and saves the results to the disk

% Create an empty error message
    ME20=[];

% Try to run the function, otherwise print the error message
try

    % Convert the filepath from string to char
        inputFilePath = convertStringsToChars(inputFilePath);

    % Load the relevant variables
        load(inputFilePath,'cellType_OverTime')
        load(inputFilePath,'timeInMin')
        load(inputFilePath,'micronsPerPixel_X')

% Perform the analysis on cell type 1
    % Create a structure to store the data
        domainWidthinMicrons_RepByPixel_CellType1 = struct;
    % Choose the timepoints to analyze
        % All timepoints
            domainWidthinMicrons_RepByPixel_CellType1.tp2Analze = 1:size(cellType_OverTime,3);
        % Every 100th timepoint
            %domainWidthinMicrons_RepByPixel.tp2Analze = 1:100:size(cellType_onGrid_OverTime,3);
        % Every 30 minutes
            % % Calculate the target times
            %     timeInMin_Target = 0:30:max(timeInMin);
            % % Calculate the time different between all possible pairs of
            % % the target time and the actual timepoints
            %     timeDifferenceAllPairs = abs(timeInMin(:)-timeInMin_Target(:).');
            % % Find the closest simulated timepoint to each experimental
            % % timepoint
            %     [M, I] = min(timeDifferenceAllPairs,[],1);
            % % Keep only the unique timepoints
            %     domainWidthinMicrons_RepByPixel.tp2Analze = unique(I');
    % Choose the cell diameter
        domainWidthinMicrons_RepByPixel_CellType1.cellDiameterInUM = 16;    
    % Input the target spatial resolution
        domainWidthinMicrons_RepByPixel_CellType1.micronsPerPixel_X_Target = micronsPerPixel_X;
    % Input the current spatial resolution
        domainWidthinMicrons_RepByPixel_CellType1.micronsPerPixel_X_Original = micronsPerPixel_X;
    % Choose which cell type to analyze
        cellTypeID = 1;
    % Choose whether to plot the data
        doPlot = true;

    % Run the fitting algorithm
        [domainWidthinMicrons_RepByPixel_CellType1.Mean, domainWidthinMicrons_RepByPixel_CellType1.STD, domainWidthinMicrons_RepByPixel_CellType1.SE, ...
            domainWidthinMicrons_RepByPixel_CellType1.Skew, domainWidthinMicrons_RepByPixel_CellType1.Kurtosis] = ...
                calculateDomainSize_SkelDist_SizeDeptFill(cellType_OverTime(:,:,domainWidthinMicrons_RepByPixel_CellType1.tp2Analze), ...
                domainWidthinMicrons_RepByPixel_CellType1.cellDiameterInUM, domainWidthinMicrons_RepByPixel_CellType1.micronsPerPixel_X_Original, ...
                domainWidthinMicrons_RepByPixel_CellType1.micronsPerPixel_X_Target,cellTypeID, doPlot);

    % Save the analyzed data to the output file
        save(inputFilePath,'domainWidthinMicrons_RepByPixel_CellType1','-append');


% Perform the analysis on cell type 2
    % Create a structure to store the data
        domainWidthinMicrons_RepByPixel_CellType2 = struct;
    % Choose the timepoints to analyze
        % All timepoints
            domainWidthinMicrons_RepByPixel_CellType2.tp2Analze = 1:size(cellType_OverTime,3);
        % Every 100th timepoint
            %domainWidthinMicrons_RepByPixel.tp2Analze = 1:100:size(cellType_onGrid_OverTime,3);
        % Every 30 minutes
            % % Calculate the target times
            %     timeInMin_Target = 0:30:max(timeInMin);
            % % Calculate the time different between all possible pairs of
            % % the target time and the actual timepoints
            %     timeDifferenceAllPairs = abs(timeInMin(:)-timeInMin_Target(:).');
            % % Find the closest simulated timepoint to each experimental
            % % timepoint
            %     [M, I] = min(timeDifferenceAllPairs,[],1);
            % % Keep only the unique timepoints
            %     domainWidthinMicrons_RepByPixel.tp2Analze = unique(I');
    % Choose the cell diameter
        domainWidthinMicrons_RepByPixel_CellType2.cellDiameterInUM = 16;    
    % Input the target spatial resolution
        domainWidthinMicrons_RepByPixel_CellType2.micronsPerPixel_X_Target = micronsPerPixel_X;
    % Input the current spatial resolution
        domainWidthinMicrons_RepByPixel_CellType2.micronsPerPixel_X_Original = micronsPerPixel_X;
    % Choose which cell type to analyze
        cellTypeID = 2;
    % Choose whether to plot the data
        doPlot = true;

    % Run the fitting algorithm
        [domainWidthinMicrons_RepByPixel_CellType2.Mean, domainWidthinMicrons_RepByPixel_CellType2.STD, domainWidthinMicrons_RepByPixel_CellType2.SE, ...
            domainWidthinMicrons_RepByPixel_CellType2.Skew, domainWidthinMicrons_RepByPixel_CellType2.Kurtosis] = ...
                calculateDomainSize_SkelDist_SizeDeptFill(cellType_OverTime(:,:,domainWidthinMicrons_RepByPixel_CellType2.tp2Analze), ...
                domainWidthinMicrons_RepByPixel_CellType2.cellDiameterInUM, domainWidthinMicrons_RepByPixel_CellType2.micronsPerPixel_X_Original, ...
                domainWidthinMicrons_RepByPixel_CellType2.micronsPerPixel_X_Target,cellTypeID, doPlot);

    % Save the analyzed data to the output file
        save(inputFilePath,'domainWidthinMicrons_RepByPixel_CellType2','-append');

catch errorMessage
    % Record the error message
        ME20 = errorMessage;
    % Save the error message to the output file
        save(inputFilePath,'ME20','-append');
end

end
