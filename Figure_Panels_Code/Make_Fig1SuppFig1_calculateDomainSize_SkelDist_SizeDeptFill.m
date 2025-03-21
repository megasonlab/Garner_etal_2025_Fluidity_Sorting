
%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/19

% Clear the system
    close all;
    clear all;

% To run this code, update the file location to wherever you want to store
% the output figure files
    figOutputFolderPath = ['\Garner_etal_2025_Fluidity_Sorting\Figures\'];
    expmtMasterFolderPath = ['\Garner_etal_2025_Fluidity_Sorting\Experimental_Data\'];
    
%% Choose the data to plot
    % Pull out the folder
        dataFolderPath = [expmtMasterFolderPath 'Processed_Data\240420_spheroid_assay\'];
    % Pull out the mat files
        matFiles = dir([dataFolderPath '*_MaxIP_Segmentation.mat']);
    % Sort the files the way a human would
       [matFiles] = natsortfiles(matFiles);
    % Pull out the filename
       inputFilePath = [dataFolderPath matFiles(3).name];

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
            % Calculate the target times
                timeInMin_Target = 18*60;
            % Calculate the time different between all possible pairs of
            % the target time and the actual timepoints
                timeDifferenceAllPairs = abs(timeInMin(:)-timeInMin_Target(:).');
            % Find the closest simulated timepoint to each experimental
            % timepoint
                [M, I] = min(timeDifferenceAllPairs,[],1);
            % Keep only the unique timepoints
                domainWidthinMicrons_RepByPixel_CellType1.tp2Analze = unique(I');
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
        cellType_onGrid_OverTime = cellType_OverTime(:,:,domainWidthinMicrons_RepByPixel_CellType1.tp2Analze);
        cellDiameterInUM = domainWidthinMicrons_RepByPixel_CellType1.cellDiameterInUM;
        micronsPerPixel_X_Original = domainWidthinMicrons_RepByPixel_CellType1.micronsPerPixel_X_Original;
        micronsPerPixel_X_Target = domainWidthinMicrons_RepByPixel_CellType1.micronsPerPixel_X_Target;

% Prepare for image upsampling
    if micronsPerPixel_X_Target~=micronsPerPixel_X_Original
        % Calculate the number pixels in the simulated image when equivalently
        % scaled to the experimental spatial resolution 
            newNumPixelsPerSide = round(size(cellType_onGrid_OverTime,1)*micronsPerPixel_X_Original/micronsPerPixel_X_Target);
        % Update the microns per pixel
            micronsPerPixel_X = micronsPerPixel_X_Original*size(cellType_onGrid_OverTime,1)/newNumPixelsPerSide;
    else
        micronsPerPixel_X = micronsPerPixel_X_Original;
    end

% Calculate the number of pixels per cell
    numPixelsPerCell = round(pi*(cellDiameterInUM/2)^2/(micronsPerPixel_X^2));
    widthPixelsPerCell = round(cellDiameterInUM/micronsPerPixel_X);

% Pre-allocate space to store the results
    domainWidthinMicrons_RepByPixel_Mean = nan([size(cellType_onGrid_OverTime,3) 1]);
    domainWidthinMicrons_RepByPixel_STD = nan([size(cellType_onGrid_OverTime,3) 1]);
    domainWidthinMicrons_RepByPixel_SE = nan([size(cellType_onGrid_OverTime,3) 1]);
    domainWidthinMicrons_RepByPixel_Skew = nan([size(cellType_onGrid_OverTime,3) 1]);
    domainWidthinMicrons_RepByPixel_Kurtosis = nan([size(cellType_onGrid_OverTime,3) 1]);

for tpNum = 1:size(cellType_onGrid_OverTime,3)

        % Pull out the cell type grid for this timepoint
            cellType_OnGrid_ThisTP = cellType_onGrid_OverTime(:,:,tpNum);

        % Perform the width calculation
            % Pull out the binary image for a single cell type
                BW_Image1 = false(size(cellType_OnGrid_ThisTP));
                BW_Image1(cellType_OnGrid_ThisTP==cellTypeID) = true;
            % Resize the image 
            if micronsPerPixel_X_Target~=micronsPerPixel_X_Original
                BW_Image1 = imresize(BW_Image1,[newNumPixelsPerSide newNumPixelsPerSide],'nearest');
            end
        % Find the initial domain size
            % Fill small holes
                BW_Image2 = ~bwareaopen(~BW_Image1, numPixelsPerCell,4);
            % Remove objects smaller than 1 cell
                BW_Image2 = bwareafilt(BW_Image2,[0.75*numPixelsPerCell Inf]);
            % Smooth object edges by ~ 1 cell width
                windowSizeInPixels = widthPixelsPerCell;
                kernel = ones(windowSizeInPixels) / windowSizeInPixels ^ 2;
                BW_Image2 = conv2(single(BW_Image2), kernel, 'same');
                BW_Image2 = BW_Image2 > 0.4; 
            % Perform the skeletonization
                skel = bwskel(BW_Image2);
            % Perform the distrance transform (distance of each pixel to a "false" pixel)
                D_Pix = bwdist(~BW_Image2);
                D_Pix(D_Pix~=0) = D_Pix(D_Pix~=0)-0.5; 
            % Calculate the domain widths along the skeleton
                domainWidthinMicrons_SkelDist = D_Pix(skel(:))*micronsPerPixel_X*2;
            % Weight each skeleton point by its domain size (so each pixel
            % along that width gets a fair say)
                domainWidthinMicrons_SkelDist_RepByPixel = repelem(domainWidthinMicrons_SkelDist,round(2*D_Pix(skel(:))));

        % Plot the results
            if doPlot
                % Calculate the spatial position vectors
                    xValsInMicrons = ((1:size(BW_Image1,1))-0.5)*micronsPerPixel_X;
                    yValsInMicrons = ((1:size(BW_Image1,2))-0.5)*micronsPerPixel_X;
                % Calculate the possible domain sizes
                    possibleDomainSizesInUM = (0:(0.5*widthPixelsPerCell):(sqrt(2)*size(BW_Image1,1)));
                figure(1)            
                subplot(2,3,1)
                imagesc(xValsInMicrons,yValsInMicrons,BW_Image1)
                xlabel('X-position (μm)')
                ylabel('Y-position (μm)')
                title('Original Binary Image')
                subplot(2,3,2)
                imagesc(xValsInMicrons,yValsInMicrons,BW_Image2)
                xlabel('X-position (μm)')
                ylabel('Y-position (μm)')
                title({'Small holes filled','Objects < 1 cell removed','Edges smoothed by 1 cell diameter'})
                subplot(2,3,3)
                imagesc(xValsInMicrons,yValsInMicrons,BW_Image2 + imdilate(skel,strel("square",5)))
                xlabel('X-position (μm)')
                ylabel('Y-position (μm)')
                title('Skeletonization') 
                subplot(2,3,4)
                imagesc(xValsInMicrons,yValsInMicrons,D_Pix.*micronsPerPixel_X)
                colorbar
                xlabel('X-position (μm)')
                ylabel('Y-position (μm)')
                title('Distance to edge (μm)') 
                subplot(2,3,5)
                imagesc(xValsInMicrons,yValsInMicrons,imdilate(D_Pix.*skel.*micronsPerPixel_X,strel("square",5)),...
                    'alphadata',imdilate(D_Pix.*skel.*micronsPerPixel_X,strel("square",5))>0)
                colorbar
                xlabel('X-position (μm)')
                ylabel('Y-position (μm)')
                title('Distance of skeleton to edge (μm)')  
                subplot(2,3,6)
                histogram(domainWidthinMicrons_SkelDist_RepByPixel,...
                    possibleDomainSizesInUM,'Normalization','pdf')
                legend(sprintf('Domain size \n %1.1f +- %1.1f μm',...
                    mean(domainWidthinMicrons_SkelDist_RepByPixel),...
                    std(domainWidthinMicrons_SkelDist_RepByPixel)))
                xlim([0 ceil(max(domainWidthinMicrons_SkelDist_RepByPixel)/100)*100])
                title('Distribution of domain sizes (μm)')
                xlabel('Domain size (μm)')
                ylabel('PDF')
                drawnow;
            end
            

        % Loop through each object, fill holes smaller than the domain
        % size, and redo smoothing based on the mean domain size
            % Fill small holes
                BW_Image2 = ~bwareaopen(~BW_Image1, numPixelsPerCell,4);
            % Remove objects smaller than 1 cell
                BW_Image2 = bwareafilt(BW_Image2,[0.75*numPixelsPerCell Inf]);
            % Label and count the objects
                [labeledImage, numObjects] = bwlabel(BW_Image2);
            % Plot the image
            % if doPlot
            %     figure(2)
            %     subplot(1,3,1)
            %     imagesc(labeledImage)
            % end
            % Create a new BW image
                BW_Image3 = false(size(BW_Image1));
            for objNum = 1:numObjects
                % Extract the BW object
                    BW_Obj = (labeledImage==objNum);
                % Extract the object's skeleton
                    skelObj = (skel&BW_Obj);
                if any(skelObj(:))
                    % % Plot the data
                    % if doPlot
                    %     subplot(1,3,2)
                    %     imagesc(imdilate(D_Pix.*skelObj,strel('disk',5)))
                    %     colorbar                
                    %     title({'Shortest distance to','the nearest black pixel (pixels)'})  
                    % end
                    % Calculate the mean distance to the edge in pixels
                        meanDomainSizeObjPix = double(round(mean(D_Pix(skelObj))));
                    % Fill holes smaller than the mean distance^2
                        BW_Obj2 = ~bwareaopen(~BW_Obj, round(meanDomainSizeObjPix.^2), 4);
                    % Smooth object edges by the mean distance (or the cell
                    % size, whichever is smaller)
                        windowSizeInPixels = max([meanDomainSizeObjPix,widthPixelsPerCell]);
                        kernel = ones(windowSizeInPixels) / windowSizeInPixels ^ 2;
                        BW_Obj2 = conv2(single(BW_Obj2), kernel, 'same');
                        BW_Obj2 = BW_Obj2 > 0.4; 
                    % Update the BW object
                        BW_Image3(BW_Obj2) = true;
                else
                    BW_Image3(BW_Obj) = true;
                end
                % % Plot the result
                % if doPlot
                %     subplot(1,3,3)
                %     imagesc(bwlabel(BW_Image3))
                % end
            end      
          
        % Perform the width calculation again
            % Perform the skeletonization
                skel = bwskel(BW_Image3);
            % Perform the distrance transform (distance of each pixel to a "false" pixel)
                D_Pix = bwdist(~BW_Image3);
                D_Pix(D_Pix~=0) = D_Pix(D_Pix~=0)-0.5; 

        % Separate the solid+square objects (for which skeletonization
        % fails) to calculate the domain size separately for these pixels
            % Classify the objects
                % Label and count the objects
                    [labeledImage, numObjects] = bwlabel(BW_Image3);
                % Calculate the BW object properties
                    BW_Image3_props = regionprops(labeledImage,'MajorAxisLength','MinorAxisLength','Solidity','Area');
                    majorAxisLength = [BW_Image3_props.MajorAxisLength];
                    minorAxisLength = [BW_Image3_props.MinorAxisLength];
                    symmetryScore = minorAxisLength./majorAxisLength;
                    solidity = [BW_Image3_props.Solidity];
                    area_Pix = [BW_Image3_props.Area];
                % Find the most solid and symmetric object (for which the
                % skeletonization fails)
                    solid_square_Obj = find((symmetryScore>0.85)&(solidity>0.9));
                % Create two separate BW images containing the solid+symmetric
                % objects
                    BW_Image3_solid_square = false(size(labeledImage));
                    BW_Image3_solid_square(ismember(labeledImage(:),solid_square_Obj)) = true;
                    BW_Image3_elongated = BW_Image3;
                    BW_Image3_elongated(ismember(labeledImage(:),solid_square_Obj)) = false;
                % Plot the objects to show that it worked
                % if doPlot
                %     figure(3)
                %     subplot(1,3,1)
                %     imagesc(labeledImage)
                %     title('All objects')
                %     subplot(1,3,2)
                %     imagesc(BW_Image3_solid_square)
                %     title('Solid symmetric objects')
                %     subplot(1,3,3)
                %     imagesc(BW_Image3_elongated)
                %     title('All other objects')
                % end
            % For the solid+square objects, use the mean of the major and 
            % minor axes lengths as the domain size, and weight by the 
            % number of pixels in the object area
                % Calculate the domain size for these objects
                    domainVals_Pix = minorAxisLength(solid_square_Obj);
                % Count the number of pixels in each object
                    repVals = area_Pix(solid_square_Obj);
                % Create the domain values
                    domainWidthinMicrons_solid_symmetric_RepByPixel = ...
                        repelem(domainVals_Pix.*micronsPerPixel_X,repVals)';
            % For the non solid+square objects, use the skeleton
                % Calculate the domain widths along the skeleton
                    domainWidthinMicrons_SkelDist = D_Pix(skel(:)&BW_Image3_elongated(:))*micronsPerPixel_X*2;
                % Weight each skeleton point by its domain size (so each pixel
                % along that width gets a fair say)
                    domainWidthinMicrons_SkelDist_RepByPixel = repelem(domainWidthinMicrons_SkelDist,round(2*D_Pix(skel(:)&BW_Image3_elongated(:))));
            % Concatenate the domain size vectors
                domainWidthinMicrons_RepByPixel = [domainWidthinMicrons_SkelDist_RepByPixel ; domainWidthinMicrons_solid_symmetric_RepByPixel];

        % Record the results
            domainWidthinMicrons_RepByPixel_Mean(tpNum) = mean(domainWidthinMicrons_RepByPixel);
            domainWidthinMicrons_RepByPixel_STD(tpNum) = std(domainWidthinMicrons_RepByPixel);
            domainWidthinMicrons_RepByPixel_SE(tpNum) = std(domainWidthinMicrons_RepByPixel)/sqrt(length(domainWidthinMicrons_RepByPixel));
            domainWidthinMicrons_RepByPixel_Skew(tpNum) = skewness(domainWidthinMicrons_RepByPixel);
            domainWidthinMicrons_RepByPixel_Kurtosis(tpNum) = kurtosis(domainWidthinMicrons_RepByPixel);

        if doPlot
        % Plot the results
            figure(4)            
            subplot(2,3,1)
            imagesc(xValsInMicrons,yValsInMicrons,BW_Image1)
            xlabel('X-position (μm)')
            ylabel('Y-position (μm)')
            title('Original Binary Image')
            set(gca,'FontName','Helvetica','FontSize',5)
            box on;
            subplot(2,3,2)
            imagesc(xValsInMicrons,yValsInMicrons,BW_Image3)
            xlabel('X-position (μm)')
            ylabel('Y-position (μm)')
            title({'Small holes filled','Objects < 1 cell removed','Edges smoothed by 1 cell diameter'})
            set(gca,'FontName','Helvetica','FontSize',5)
            box on;
            subplot(2,3,3)
            imagesc(xValsInMicrons,yValsInMicrons,BW_Image3 + imdilate(skel.*BW_Image3_elongated,strel("square",5)))
            xlabel('X-position (μm)')
            ylabel('Y-position (μm)')
            title('Skeletonization') 
            set(gca,'FontName','Helvetica','FontSize',5)
            box on;
            subplot(2,3,4)
            imagesc(xValsInMicrons,yValsInMicrons,D_Pix.*micronsPerPixel_X)
            colorbar
            xlabel('X-position (μm)')
            ylabel('Y-position (μm)')
            title('Distance to edge (μm)') 
            set(gca,'FontName','Helvetica','FontSize',5)
            box on;
            axis equal;
            subplot(2,3,5)
            imagesc(xValsInMicrons,yValsInMicrons,imdilate(D_Pix.*skel.*BW_Image3_elongated.*micronsPerPixel_X,strel("square",5)),...
                'alphadata',imdilate(D_Pix.*skel.*BW_Image3_elongated.*micronsPerPixel_X,strel("square",5))>0)
            colorbar
            xlabel('X-position (μm)')
            ylabel('Y-position (μm)')
            title('Distance of skeleton to edge (μm)')  
            set(gca,'FontName','Helvetica','FontSize',5)
            box on;
            axis equal;
            subplot(2,3,6)
            histogram(domainWidthinMicrons_RepByPixel,...
                possibleDomainSizesInUM,'Normalization','pdf')
            legend(sprintf('Domain size \n %1.1f +- %1.1f μm',...
                mean(domainWidthinMicrons_RepByPixel),...
                std(domainWidthinMicrons_RepByPixel)))
            xlim([0 ceil(max(domainWidthinMicrons_RepByPixel)/100)*100])
            title('Distribution of domain sizes (μm)')
            xlabel('Domain size (μm)')
            ylabel('PDF')
            drawnow;
            set(gca,'FontName','Helvetica','FontSize',5)
            box on;
        end

end


     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig1SuppFig1_DomainSizeCalc'];
        % Choose the figure size and position
            figHandle = figure(4);
            figHandle.Position =  [250   375   500   300];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])