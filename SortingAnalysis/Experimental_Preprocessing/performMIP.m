%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This script performs a maximum intensity projection (MIP) along the 
% z-axis of all image volumes included in the specified nd2 file, 
% saving the MIP as a separate tiff for each series.

%% Clear the system
    close all;
    clear all;

%% User selections

% Input the image file path
    expmtMasterFolderPath = ['\Garner_etal_2025_Fluidity_Sorting\Experimental_Data\'];
    imageFilePath = [expmtMasterFolderPath 'Raw_Data\240420_spheroid_assay\L929_2h_220um_timeseries.nd2']; 
    %imageFilePath = [expmtMasterFolderPath 'Raw_Data\240418_spheroid_assay\L929_2h_220um_timeseries.nd2']; 

% Choose whether to delete old MIP files and redo the MIP 
    redoAnalysis = true;

%% Load the image information

% Use bioformats tools to read the file (but do not load the images)
    reader = bfGetReader(imageFilePath);

% Use bioformats tools to pull out the relevant metadata
    
    % Read the metadata
        omeMeta = reader.getMetadataStore();

    % Get the data type
        dtype = omeMeta.getPixelsType(0).getValue();
    
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
        % Number of channels
            numChannels = omeMeta.getPixelsSizeC(0).getValue();
        % Number of series
            numSeries = omeMeta.getImageCount();
    
    % Get the pixel sizes in microns
        % Pull the microns per pixel from java and convert to double
            micronsPerPixel_X = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();
            micronsPerPixel_Y = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();
            micronsPerPixel_Z = omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROMETER).doubleValue();

%% Loop through each series and timepoint and calculate the MIP

for iSeries = 1:numSeries

    % Print the series number
        fprintf('Analyzing series %d of %d.\n\n',iSeries,numSeries)

    % Choose where to save the maximum intensity projection
        MIPFilePath = replace(imageFilePath,'Raw_Data','Processed_Data');
        MIPFilePath = [MIPFilePath(1:(find(MIPFilePath=='.',1,'last'))-1) ...
            sprintf('_S%1i_',iSeries) 'MaxIP_MATLAB.tiff'];

    % Delete the file if it already exists and the user has chosen to
    % re-do the analysis
    if (exist(MIPFilePath, 'file')==2)&(redoAnalysis)
        delete(MIPFilePath);
    end
        
    % If the file does not already exist, run the analysis
    if exist(MIPFilePath, 'file')~=2
            
        %  Use bioformats tools to select the appropriate series
            reader.setSeries(iSeries - 1);

        % Preallocate space to store the image MIPs
            imageVol_MIP = zeros(numPixels_X,numPixels_Y,1,numChannels,numTimePoints,char(dtype));
    
        for iT = 1:numTimePoints
        
            % Print the timepoint
                fprintf('Analyzing timepoint %d of %d.\n',iT,numTimePoints)
        
            % Load the image volume
                % Loop through each channel and load the image
                for iC = 1:numChannels
                    % Preallocate space to store the image volume
                        imageVol = nan(numPixels_X,numPixels_Y,numPixels_Z);                    
                    % Loop through each z-slice and load the image
                    for iZ = 1:numPixels_Z
                            % Use bioformats tools to pull out the appropriate plane
                                iPlane = reader.getIndex(iZ - 1, iC - 1, iT - 1) + 1;
                            % Use bioformats tools to load the z-plane 
                            % and save it to the image volume
                                imageVol(:,:,iZ) = bfGetPlane(reader, iPlane);
                    end
                    % Perform the max intensity projection along z
                        imageVol_MIP(:,:,1,iC,iT) = squeeze(max(imageVol,[],3));   
                end
        
        end
  
        % Use bioformats tools to save the MIP
            bfsave(imageVol_MIP,MIPFilePath, 'dimensionOrder', 'XYZCT')

    end

end
