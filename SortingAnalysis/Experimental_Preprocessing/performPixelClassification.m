%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This script performs a pixel-based segmentation of the MIP images into
% cell type #1, cell type #2, and background 

% Clear the system
    close all;
    clear all;

% To run tihs code, update the following path to the data on your end   
    expmtMasterFolderPath = ['\Garner_etal_2025_Fluidity_Sorting\Experimental_Data\'];

% Choose whether to delete old analysis files and redo the pixel classification 
    redoAnalysis = false;

%% User selections

% 2024/04/20
    % Input the image file path
        imageFilePath = [expmtMasterFolderPath 'Raw_Data\240420_spheroid_assay\L929_2h_220um_timeseries.nd2']; 
    % Choose the channels to analyze
        cellType1ChannelIdx = 1;
        cellType2ChannelIdx = 2;
    % Choose the thresholds for each series 
        thresholdLevels = [0.002, 0.0035, 0.002, 0.0035, 0.0035, 0.002, 0.0035, 0.002];
    % Input the time delay between plating and imaging
        timeDelayInMin = 120;
    % Input the series numbers
        seriesNums = 1:8;
    % Input the number of cells for each condition
        numCells = 80000*ones(8,1);
    % Input the number of cells for each condition
        cellType1Name = {'L929_Cdh2-GFP_low','L929_Cdh2-GFP_high','L929_Cdh3-GFP_low','L929_Cdh3-GFP_high',...
            'L929_Cdh3-GFP_high','L929_Cdh3-GFP_low','L929_Cdh2-GFP_high','L929_Cdh2-GFP_low'};
        cellType2Name = {'L929_Cdh1-RFP_low','L929_Cdh1-RFP_high','L929_Cdh1-RFP_low','L929_Cdh1-RFP_high',...
            'L929_Cdh1-RFP_high','L929_Cdh1-RFP_low','L929_Cdh1-RFP_high','L929_Cdh1-RFP_low'};
    % Input the date
        dateVal = '20240420'; 

% 2024/04/18
    % Input the image file path
        imageFilePath = [expmtMasterFolderPath 'Raw_Data\240418_spheroid_assay\L929_2h_220um_timeseries.nd2']; 
    % Choose the channels to analyze
        cellType1ChannelIdx = 1;
        cellType2ChannelIdx = 2;
    % Choose the thresholds for each series 
        thresholdLevels = [0.002, 0.0035, 0.002, 0.0035, 0.0035, 0.002, 0.0035, 0.002];
    % Input the time delay between plating and imaging
        timeDelayInMin = 120;
    % Input the series numbers
        seriesNums = 1:8;
    % Input the number of cells for each condition
        numCells = 80000*ones(8,1);
    % Input the number of cells for each condition
        cellType1Name = {'L929_Cdh2-GFP_low','L929_Cdh2-GFP_high','L929_Cdh3-GFP_low','L929_Cdh3-GFP_high',...
            'L929_Cdh3-GFP_high','L929_Cdh3-GFP_low','L929_Cdh2-GFP_high','L929_Cdh2-GFP_low'};
        cellType2Name = {'L929_Cdh1-RFP_low','L929_Cdh1-RFP_high','L929_Cdh1-RFP_low','L929_Cdh1-RFP_high',...
            'L929_Cdh1-RFP_high','L929_Cdh1-RFP_low','L929_Cdh1-RFP_high','L929_Cdh1-RFP_low'};
    % Input the date
        dateVal = '20240418'; 


%% Load the image information

    % Use bioformats tools to read the file (but do not load the images)
        reader = bfGetReader(imageFilePath);
    
    % Use bioformats tools to pull out the relevant metadata
        
        % Read the metadata
            omeMeta = reader.getMetadataStore();
        
        % Read in the dimension sizes
            % Image width in pixels
                numPixels_X = omeMeta.getPixelsSizeX(0).getValue(); 
            % Image height in pixels
                numPixels_Y = omeMeta.getPixelsSizeY(0).getValue();
            % Number of z-slices
                numPixels_Z = omeMeta.getPixelsSizeZ(0).getValue();
            % Concatenate into single array
                numPixels = [numPixels_X numPixels_Y numPixels_Z];
            % Number of time points
                numTimePoints = omeMeta.getPixelsSizeT(0).getValue();
               % numTimePoints = 49;
            % Number of channels
                numChannels = omeMeta.getPixelsSizeC(0).getValue();
            % Number of series
                numSeries = omeMeta.getImageCount();
        
        % Get the pixel sizes in microns
            % Pull the microns per pixel from java and convert to double
            % X
                micronsPerPixel_X = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();
            % Y
                micronsPerPixel_Y = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();
            % Z
                micronsPerPixel_Z = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER).doubleValue();
            % Concatenate into single array
                micronsPerPixel = [micronsPerPixel_X micronsPerPixel_Y micronsPerPixel_Z];
    
        % Get the time stamps
            % Pull out the actual time stamps
                % Pre-allocate space to store the data
                    timeStampInSec = nan(numTimePoints,1);
                % Choose the plane
                    iZ = 1;
                    iC = 1;
                    iSeries = 1;
                    reader.setSeries(iSeries - 1);
                % Loop through each timepoint and pull out the time
                try
                    for iT = 1:numTimePoints
                        % Choose the plane
                            iPlane = reader.getIndex(iZ - 1, iC - 1, iT - 1) + 1;
                        % Pull out the time stamp
                            timeStampInSec(iT) = omeMeta.getPlaneDeltaT(0,iPlane).value().doubleValue();
                    end
                catch
                    numTimePoints=iT-1;
                end

            % Convert to time step increments
                timeStampStepInSec = diff(timeStampInSec);
            % Calculate an average time step
                timeStepInSeconds = nanmean(timeStampStepInSec);

%% Loop through each series and timepoint and analyze them

% Find all the MIP files in the same folder as the raw data
    % Pull out the folder path to the MIP images
        folderPath = imageFilePath(1:(find(imageFilePath=='\',1,'last')));
        folderPath = replace(folderPath,'Raw_Data','Processed_Data');
    % Pull out the filepaths to all the MIP files in the folder
        MIPFileNames = dir([folderPath '*_MaxIP_MATLAB.tiff']);

for iSeries = 1%:numSeries

    % Print the series number
        fprintf('Analyzing series %d.\n',iSeries)

    % Update the output file paths
        % Pull out the MIP file path for this series
            MIPFilePath = [folderPath MIPFileNames(iSeries).name];
        % Choose the file path to save the figures for the domain size fit over time
            pixelClassificationFilePath = [MIPFilePath(1:(strfind(MIPFilePath,'MaxIP_MATLAB')-1)) 'MaxIP_Segmentation.tiff'];
        % Choose the file path to save the pixel classification data
            matFilePath = [pixelClassificationFilePath(1:(find(pixelClassificationFilePath=='.',1,'last'))-1) '.mat'];

    % Read the file (but do not load the images)
        % Read the file (but do not load the images)
            reader_MIP = bfGetReader(MIPFilePath);
        % Read the metadata
            omeMeta_MIP = reader_MIP.getMetadataStore();
        % Number of time points
            numTimePoints_MIP = omeMeta_MIP.getPixelsSizeT(0).getValue();

    % If the file does not already exist, run the analysis
    if (exist(pixelClassificationFilePath, 'file')==2)&(redoAnalysis)
        delete(pixelClassificationFilePath)
    end    
    if (exist(matFilePath, 'file')==2)&(redoAnalysis)
        delete(matFilePath)
    end

    if exist(pixelClassificationFilePath, 'file')~=2

        % Pre-allocate space to store the cell type matrix
            cellType_OverTime = nan(numPixels_X,numPixels_Y,numTimePoints_MIP);
    
        % Make a time vector
            timeInMin = (0:(numTimePoints_MIP-1))*timeStepInSeconds/60; 
        
        for iT = 1:numTimePoints_MIP
        
            % Print the timepoint
                fprintf('Analyzing timepoint %d of %d.\n',iT,numTimePoints_MIP)
        
            % Load the image volume
                % Pull out the appropriate planes
                    iPlane1 = reader_MIP.getIndex(iZ - 1, cellType1ChannelIdx -1, iT - 1) + 1;
                    iPlane2 = reader_MIP.getIndex(iZ - 1, cellType2ChannelIdx -1, iT - 1) + 1;
                % Load the planes and convert to double precision values from 0 to 1
                    % where 1 maps to the largest possible for that specific
                    % data type (e.g., for unit16 1=2^16-1, for unit8 1=2^8-1)
                    cellType1ChannelIdxImage_MIP = im2double(bfGetPlane(reader_MIP, iPlane1));
                    cellType2ChannelIdxImage_MIP = im2double(bfGetPlane(reader_MIP, iPlane2));
    
            % Match the histograms of the image so intensity values are
            % equivalent between the channels
            if mean(cellType1ChannelIdxImage_MIP(:))>mean(cellType2ChannelIdxImage_MIP(:))
                cellType2ChannelIdxImage_MIP = ...
                    imhistmatch(cellType2ChannelIdxImage_MIP,cellType1ChannelIdxImage_MIP,10000,'method','polynomial');
            else
                cellType1ChannelIdxImage_MIP = ...
                    imhistmatch(cellType1ChannelIdxImage_MIP,cellType2ChannelIdxImage_MIP,10000,'method','polynomial');
            end
    
            %% Segment each channel
                % Pull out the user-defined thresholds
                    level1 = thresholdLevels(iSeries);
                % Threshold the images
                    BW1 = imbinarize(cellType1ChannelIdxImage_MIP,level1);
                    BW2 = imbinarize(cellType2ChannelIdxImage_MIP,level1);
                % Fill the holes smaller than 100 pixels
                    BW1 = ~bwareaopen(~BW1, 20^2);
                    BW2 = ~bwareaopen(~BW2, 20^2);
    
            % Perform pixel-based cell type classification
                % Pre-allocate the image
                    cellType = nan(size(cellType1ChannelIdxImage_MIP));
                % For pixels that are only identified as one cell type (but not
                % both), add them to the matrix
                    cellType(BW1&(~BW2)) = 1;
                    cellType(BW2&(~BW1)) = 2;
                % For pixels identified as both cell types, choose the one with the
                % highest intensity
                    both = nan([sum(BW2(:)&BW1(:)),1]);
                    both(cellType1ChannelIdxImage_MIP(BW2&BW1)>...
                        cellType2ChannelIdxImage_MIP(BW2&BW1)) = 1;
                    both(cellType1ChannelIdxImage_MIP(BW2&BW1)<...
                        cellType2ChannelIdxImage_MIP(BW2&BW1)) = 2;
                    cellType(BW2&BW1)=both;
                % For all other pixels, set them as the mean value
                    cellType(isnan(cellType)) = 1.5; 
                % Save this in the cell idendity matrix for all timepoints
                    cellType_OverTime(:,:,iT) = cellType;
                
            % Plot the results
                figure(1)
                % Plot the MIP (intensity adjusted)
                    subplot(2,4,1)
                    imagesc((1:numPixels_X)*micronsPerPixel_X,...
                        (1:numPixels_X)*micronsPerPixel_X,...
                        imadjust(cellType1ChannelIdxImage_MIP))
                    ax1 = gca;
                    xlabel('X-position (μm)')
                    ylabel({cellType1Name{iSeries},'','Y-position (μm)'}, 'Interpreter', 'none')
                    title('Contrast-adjusted image')
                    colorbar            
                    subplot(2,4,5)
                    imagesc((1:numPixels_X)*micronsPerPixel_X,...
                        (1:numPixels_X)*micronsPerPixel_X,...
                        imadjust(cellType2ChannelIdxImage_MIP))
                    ax2 = gca;
                    xlabel('X-position (μm)')
                    ylabel({cellType2Name{iSeries},'','Y-position (μm)'}, 'Interpreter', 'none')
                    title('Contrast-adjusted image')
                    colorbar
                % Plot the binary segmentation
                    subplot(2,4,2)
                    imagesc((1:numPixels_X)*micronsPerPixel_X,...
                        (1:numPixels_X)*micronsPerPixel_X,...
                        BW1)
                    ax3 = gca;
                    xlabel('X-position (μm)')
                    ylabel('Y-position (μm)')
                    title('Thresholded image')
                    colorbar
                    subplot(2,4,6)
                    imagesc((1:numPixels_X)*micronsPerPixel_X,...
                        (1:numPixels_X)*micronsPerPixel_X,...
                        BW2);
                    ax4 = gca;
                    xlabel('X-position (μm)')
                    ylabel('Y-position (μm)')
                    title('Thresholded image')
                    colorbar
                % Plot the histogram and threshold
                    subplot(2,4,3)
                    histogram(cellType1ChannelIdxImage_MIP(:),1000)
                    yl = ylim();
                    hold on;
                    plot([level1 level1],yl,'r-')
                    hold off;
                    xlabel('Pixel intensity')
                    ylabel('Number of pixels')
                    title('Pixel intensity distribution')
                    legend({'Pixels','Threshold'})
                    subplot(2,4,7)
                    histogram(cellType2ChannelIdxImage_MIP(:),1000)
                    yl = ylim();
                    hold on;
                    plot([level1 level1],yl,'r-')
                    hold off;
                    xlabel('Pixel intensity')
                    ylabel('Number of pixels')
                    title('Pixel intensity distribution')
                    legend({'Pixels','Threshold'})
                % Plot the segmentation    
                    subplot(2,4,4)
                    imagesc((1:numPixels_X)*micronsPerPixel_X,...
                        (1:numPixels_X)*micronsPerPixel_X,...
                        cellType_OverTime(:,:,iT))
                    ax7 = gca;
                    colorbar            
                    title('Cell Type Classification')
                    xlabel('X-position (μm)')
                    ylabel('Y-position (μm)')
                    colormap(gca,[0 1 0; 0 0 0; 1 0 1])
                % Plot the rosegmentation    
                    % Create an associated RGB image
                        MIP_RGB = nan(numPixels_X,numPixels_Y,3);
                        MIP_RGB(:,:,2) = imadjust(cellType1ChannelIdxImage_MIP);
                        MIP_RGB(:,:,1) = imadjust(cellType2ChannelIdxImage_MIP);
                        MIP_RGB(:,:,3) = imadjust(cellType2ChannelIdxImage_MIP);
                    subplot(2,4,8)
                    imagesc((1:numPixels_X)*micronsPerPixel_X,...
                        (1:numPixels_X)*micronsPerPixel_X,...
                        MIP_RGB)
                    ax8 = gca;
                    colorbar            
                    title('Original image')
                    xlabel('X-position (μm)')
                    ylabel('Y-position (μm)')
                % Link the axes
                    %linkaxes([ax1, ax2, ax3, ax4, ax7, ax8])
                % Add a title
                    sgtitle(sprintf('%s - S%i - %5i cells - %s - %s, T = %1.1f hrs',...
                        dateVal,iSeries,numCells(iSeries),cellType1Name{iSeries},...
                        cellType2Name{iSeries},timeStampInSec(iT)./3600),'Interpreter','None')
                %% Display the figure
                 %   drawnow
     
            % User has the option to stop it
            if (iT==1)
                %ginput()
            end
            
            % Save the figure to a time-lapse tiff stack       
                % Choose the figure size and position
                    figHandle = figure(1);
                    figHandle.Position =  [20   150   1500   500];   
                % Record the figure as an image
                    imagewd = getframe(gcf);            
                % Save the fitting information
                    imwrite(imagewd.cdata,pixelClassificationFilePath,'Compression','none','WriteMode', "append");
        
        end

        % Save the data to a structure
            outputVars = struct();
            outputVars.cellType_OverTime = cellType_OverTime;
            outputVars.timeInMin = timeInMin;
            outputVars.micronsPerPixel_X = micronsPerPixel_X;

        % Save the data to a new output file
            save(matFilePath,'-struct','outputVars','-v7.3');
    end

end

toc