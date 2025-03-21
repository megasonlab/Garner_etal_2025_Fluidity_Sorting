
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
    simMasterFolderPath = ['\Garner_etal_2025_Fluidity_Sorting\Simulation_Data\'];
    expmtMasterFolderPath = ['\Garner_etal_2025_Fluidity_Sorting\Experimental_Data\'];

%% Load an example dataset

% Choose the path to the data
    simDataFolderPath = [simMasterFolderPath 'Fixed_Time\Results20240408\'];

% Pull out the mat files
    matFiles = dir([simDataFolderPath '*_out.mat']);
    % Sort the files the way a human would
       [matFiles] = natsortfiles(matFiles);

% Load the number of timepoints
    load([simDataFolderPath matFiles(1).name],'numTP2Save');

% For recording the parameter inputs
    kT_All = nan([1 length(matFiles)]);  
    kT_All_Idx = nan([1 length(matFiles)]);  
    E_homo_All = nan([1 length(matFiles)]);  
    E_homo_All_Idx = nan([1 length(matFiles)]);  
    v_All = nan([1 length(matFiles)]);  
    v_All_Idx = nan([1 length(matFiles)]);  
    repNum_All = nan([1 length(matFiles)]); 

% For recording the sorting dynamics
    timeInMin_All = nan([length(matFiles) numTP2Save]);
    domainSizeInMicrons_All = nan([length(matFiles) numTP2Save]);
    domainSizeInMicronsSD_All = nan([length(matFiles) numTP2Save]);
    domainSizeInMicronsCV_All = nan([length(matFiles) numTP2Save]);

% Record how long it takes to load the data
tic

% Loop through each file and fit the domain size over time and the domain
% growth law
for matFileNum = 1:length(matFiles)

    % Print the file number
      %  matFileNum

    % Pull out the file path
        matFilePath = [simDataFolderPath matFiles(matFileNum).name];

    try

    % Load the file
        clear globalInfo
        clear parameterValsNum
        clear replicateNum
        clear timeInMin_OverTime;
        clear numTP2Save;
        clear domainWidthinMicrons_RepByPixel
        load(matFilePath,'globalInfo','parameterValsNum','replicateNum',...
            'timeInMin_OverTime','numTP2Save','domainWidthinMicrons_RepByPixel')

    % Save the parameter values
        kT_All(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));
        kT_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,3);
        E_homo_All(matFileNum) = globalInfo.parameterVals1(globalInfo.paramCombos(parameterValsNum,1));
        E_homo_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,1);
        v_All(matFileNum) = globalInfo.parameterVals2(globalInfo.paramCombos(parameterValsNum,2));
        v_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,2);
        repNum_All(matFileNum) = replicateNum;

    % Record the domain size information 
        timeInMin_All(matFileNum,:) = timeInMin_OverTime(1:numTP2Save)';
        clear domainWidthinMicrons_RepByPixel
        load(matFilePath,'domainWidthinMicrons_RepByPixel')
        domainSizeInMicrons_All(matFileNum,:) = domainWidthinMicrons_RepByPixel.Mean';
        domainSizeInMicronsSD_All(matFileNum,:) = domainWidthinMicrons_RepByPixel.STD';
        domainSizeInMicronsCV_All(matFileNum,:) = domainWidthinMicrons_RepByPixel.STD'./domainWidthinMicrons_RepByPixel.Mean';

    end


end

% Load the the size of the simulation grid in # cells
    load(matFilePath,'sizeGrid','cellRadius')

toc

%% Plot the sorting at 48 hours as a heatmap

% Calcluate the mean, SD, and CV of the domain size at 48 hrs for all
% simulation parameters
    domainSizeInMicrons_48hr = nanmean(domainSizeInMicrons_All(:,(end-50):end),2);
    domainSizeInMicronsSD_48hr = nanmean(domainSizeInMicronsSD_All(:,(end-50):end),2);
    domainSizeInMicronsCV_48hr = nanmean(domainSizeInMicronsCV_All(:,(end-50):end),2);

% Calculate the governing parameter ratios
    ERatio = kT_All;
    VRatio = E_homo_All./v_All;

% Identify the unique parameter values
    uniqueERatios = unique(ERatio);
    numUniqueERatios = length(uniqueERatios);
    uniqueVRatios = unique(VRatio);
    numUniqueVRatios = length(uniqueVRatios);

% Create a matrix to store the domain sizes
    sortState_48hr_heatmap = nan([numUniqueVRatios,numUniqueERatios]);
    sortStateSD_48hr_heatmap = nan([numUniqueVRatios,numUniqueERatios]);
    sortStateCV_48hr_heatmap = nan([numUniqueVRatios,numUniqueERatios]);

% Loop through each set of parameters and record the domain size
for ERatioNum = 1:numUniqueERatios
    for VRatioNum = 1:numUniqueVRatios
       
        % Pull out the mat files matching this combination
            matFiles2Plot = find((ERatio==uniqueERatios(ERatioNum))&(VRatio==uniqueVRatios(VRatioNum)));

        % Fill in the heatmap
        if any(~isnan(domainSizeInMicrons_48hr(matFiles2Plot)))
            sortState_48hr_heatmap(VRatioNum,ERatioNum) = mean(domainSizeInMicrons_48hr(matFiles2Plot),'omitnan');
            sortStateSD_48hr_heatmap(VRatioNum,ERatioNum) = mean(domainSizeInMicronsSD_48hr(matFiles2Plot),'omitnan');
            sortStateCV_48hr_heatmap(VRatioNum,ERatioNum) = mean(domainSizeInMicronsCV_48hr(matFiles2Plot),'omitnan');
        end
    end
end

% Plot the mean domain size at 48 hours as a function of the parameters 
    figure(1)
    subplot(2,4,1)
    imagesc(uniqueERatios,uniqueVRatios,log10(sortState_48hr_heatmap))
    title('Domain size at 48 hrs')
    c = colorbar;
    c.Label.String = 'Domain size (um)';
    c.TicksMode='manual';
    cTicks = log10(2.^(0:0.5:9));
    c.Ticks = cTicks;
    cTickLabels = arrayfun(@num2str, int16(10.^cTicks),'UniformOutput', 0);
    c.TickLabels = cTickLabels;
    set(gca,'YScale','log')
    xlabel(' kT_{eff}/E_{homo} (unitless)')
    ylabel(' E_{homo}/viscosity')
    set(gca,'YDir','normal') 
    set(gca,'FontName','Helvetica','FontSize',5); 


% Plot the mean domain size at 48 hours as a function of the parameters 
    figure(1)
    subplot(2,4,2)
    imagesc(uniqueERatios,uniqueVRatios,log10(sortStateSD_48hr_heatmap))
    title('Domain size SD at 48 hrs')
    c = colorbar;
    c.Label.String = 'Domain size SD (um)';
    c.TicksMode='manual';
    cTicks = log10(2.^(0:0.5:9));
    c.Ticks = cTicks;
    cTickLabels = arrayfun(@num2str, int16(10.^cTicks),'UniformOutput', 0);
    c.TickLabels = cTickLabels;
    set(gca,'YScale','log')
    xlabel(' kT_{eff}/E_{homo} (unitless)')
    ylabel(' E_{homo}/viscosity')
    set(gca,'YDir','normal') 
    set(gca,'FontName','Helvetica','FontSize',5); 

% Plot the CV of the domain size at 48 hours as a function of the parameters 
    figure(1)
    subplot(2,4,3)
    imagesc(uniqueERatios,uniqueVRatios,sortStateCV_48hr_heatmap)
    title('Domain size CV at 48 hrs')
    c = colorbar;
    c.Label.String = 'Domain size CV (um)';
    set(gca,'YScale','log')
    xlabel(' kT_{eff}/E_{homo} (unitless)')
    ylabel(' E_{homo}/viscosity')
    set(gca,'YDir','normal') 
    set(gca,'FontName','Helvetica','FontSize',5); 

% Calculate the fluidity as a function of the parameters
    % Load the relevant parameters to calculate the fluidity
        % Load the data
            clear cellRadius
            clear kT
            load(matFilePath,'cellRadius','kT')   
        % Pull out the thermal energy of the lab temperature in pN nm
            kT_L_pNnm = kT;
            L_nm = cellRadius*2;
        % Choose the number of homotypic neighbors for calculating the fluidity
            numSameCellTypeNeighbors2PlotFluidity = 2;        

    % Create a mesh of the parameters
        [uniqueERatios_mesh,uniqueVRatios_mesh] = meshgrid(uniqueERatios,uniqueVRatios);
        
    % Calculate the fluidity in min^-1
        fluidity_min = 60*10^3*((uniqueERatios_mesh.*uniqueVRatios_mesh.*kT_L_pNnm)./L_nm^2).*...
            exp(-2*numSameCellTypeNeighbors2PlotFluidity./uniqueERatios_mesh);
    
    % Pull out the values with unrealistic fluidity
        alpha = (fluidity_min<1);
        alpha = ones(size(fluidity_min));

% Plot the fluidity for each set of parameters
    subplot(2,4,4)
    imagesc(uniqueERatios,uniqueVRatios,log10(fluidity_min),'alphadata',alpha)
    title('Fluidity (50% sorted)')
    c = colorbar;
    c.Label.String = 'log_{10}(Fluidity, min^{-1})';
    colormap(gca,jet*0.9)
    % c.TicksMode='manual';
    % cTicks = log10(2.^(0:0.5:9));
    % c.Ticks = cTicks;
    % cTickLabels = arrayfun(@num2str, int16(10.^cTicks),'UniformOutput', 0);
    % c.TickLabels = cTickLabels;
    set(gca,'YScale','log')
    xlabel('kT_{eff}/E_{homo} (unitless)')
    ylabel('E_{homo}/viscosity')
    set(gca,'YDir','normal') 

%% For each of the datasets, plot the timecouse for sorting

% Choose the last timepoint in min to perform the comparison
    maxTime2FitInMin = 20;

% Create a space to store the best fit values
    ERatio_BestFitExp = nan(16,1);
    VRatio_BestFitExp = nan(16,1);
    CadCombo_BestFitExp = nan(16,1);
    Day_BestFitExp = nan(16,1);
    ReplicateWithinDay_BestFitExp = nan(16,1);
    ExpressionLevel_BestFitExp = nan(16,1);
    domainSizeInMicrons_BestFitExp = cell(16,1);
    domainSizeInMicronsSD_BestFitExp = cell(16,1);
    timeInMin_BestFitExp = cell(16,1);

% Create a counter
    seriesCounter = 0;

% Loop through each dataset
for datasetNum = 1:2

    % Select the data to fit
    if datasetNum==1
        % Input the image file path
            imageFilePath = [expmtMasterFolderPath 'Processed_Data\240418_spheroid_assay\L929_2h_220um_timeseries.nd2'];
        % Input the time delay between plating and imaging
            timeDelayInMin = -120;
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
        % Choose the dataset names
            dataNames = strcat(repmat({dateVal},length(seriesNums),1),...
                repmat({'_S'},length(seriesNums),1),...
                arrayfun(@num2str, seriesNums', 'UniformOutput', 0),...
                repmat({'_'},length(seriesNums),1),...
                arrayfun(@num2str, numCells(seriesNums)', 'UniformOutput', 0)',...
                repmat({'_'},length(seriesNums),1),...
                cellType1Name(seriesNums)',...
                repmat({'_'},length(seriesNums),1),...
                cellType2Name(seriesNums)');
    elseif datasetNum==2
        % Input the image file path
            imageFilePath = [expmtMasterFolderPath 'Processed_Data\240420_spheroid_assay\L929_2h_220um_timeseries.nd2'];
        % Input the time delay between plating and imaging
            timeDelayInMin = -120;
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
        % Choose the dataset names
            dataNames = strcat(repmat({dateVal},length(seriesNums),1),...
                repmat({'_S'},length(seriesNums),1),...
                arrayfun(@num2str, seriesNums', 'UniformOutput', 0),...
                repmat({'_'},length(seriesNums),1),...
                arrayfun(@num2str, numCells(seriesNums)', 'UniformOutput', 0)',...
                repmat({'_'},length(seriesNums),1),...
                cellType1Name(seriesNums)',...
                repmat({'_'},length(seriesNums),1),...
                cellType2Name(seriesNums)');
    end


    % Loop through each video within the dataset dataset

    for iSeries = seriesNums

        % Update a counter
            seriesCounter = seriesCounter+1;

        % Save the information
            CadCombo_BestFitExp(seriesCounter) = 1*contains(cellType1Name{iSeries},'Cdh3')+2*contains(cellType1Name{iSeries},'Cdh2');
            Day_BestFitExp(seriesCounter) = datasetNum;
            ReplicateWithinDay_BestFitExp(seriesCounter) = 1*(iSeries<=4)+2*(iSeries>4);
            ExpressionLevel_BestFitExp(seriesCounter) = 1*contains(cellType1Name{iSeries},'low')+2*contains(cellType1Name{iSeries},'high');
    
        % Select the path to the file
            expMatFilePath = [imageFilePath(1:(find(imageFilePath=='.',1,'last'))-1) ...
                sprintf('_S%1i_',iSeries) 'MaxIP_Segmentation.mat'];
    
        % Load the data
            clear domainWidthinMicrons_RepByPixel_CellType1;
            clear domainWidthinMicrons_RepByPixel_CellType2;
            clear timeInMin;
            clear micronsPerPixel_X;
            load(expMatFilePath,'timeInMin','micronsPerPixel_X',...
                'domainWidthinMicrons_RepByPixel_CellType1',...
                'domainWidthinMicrons_RepByPixel_CellType2')

        % Crop experimental data for comparison to the simulations
            lastTP2Fit = find(timeInMin<=maxTime2FitInMin*60,1,'last');
            timeInMin_Exp = timeInMin(1:lastTP2Fit)+timeDelayInMin;
            domainSizeInUM_Exp = domainWidthinMicrons_RepByPixel_CellType1.Mean(1:lastTP2Fit);
            domainSizeInUMSD_Exp = domainWidthinMicrons_RepByPixel_CellType1.STD(1:lastTP2Fit);
            domainSizeInUMCV_Exp = domainWidthinMicrons_RepByPixel_CellType1.STD(1:lastTP2Fit)./...
                domainWidthinMicrons_RepByPixel_CellType1.Mean(1:lastTP2Fit);
        
        % Determine whether the time vectors are the same for all simulation files
            eqVec = sum(~isnan(timeInMin_All),1)==...
                sum(bsxfun(@le,bsxfun(@minus,timeInMin_All,0.01*mode(timeInMin_All,1)),...
                mode(timeInMin_All,1)),1);

        % Pull out the closest timepoints in the simulations to the experimental measurements 
            if sum(eqVec)==length(eqVec)
            
                % Pull out the first non-nan values
                    matFileNum = find(~isnan(domainSizeInMicrons_All(:,2)),1,'first');
        
                % Pull out the simulation information
                    timeInMin_Sim = timeInMin_All(matFileNum,:);
        
                % Subset the experimental data to times available in the simulation
                    timeInMin_Exp_2Fit = timeInMin_Exp(timeInMin_Exp<=max(timeInMin_Sim));
                    domainSizeInUM_Exp_2Fit = domainSizeInUM_Exp(timeInMin_Exp<=max(timeInMin_Sim));
                    domainSizeInUMSD_Exp_2Fit = domainSizeInUMSD_Exp(timeInMin_Exp<=max(timeInMin_Sim));
                    domainSizeInUMCV_Exp_2Fit = domainSizeInUMCV_Exp(timeInMin_Exp<=max(timeInMin_Sim));
            
                % Find the associated time in the simulations that are closest to the
                % timepoints in the experiments
                    % Calculate the time different between all possible pairs of
                    % the experimental time and the actual timepoints
                        timeDifferenceAllPairs = abs(timeInMin_Sim(:)-timeInMin_Exp_2Fit(:).');
                    % Find the closest simulated timepoint to each experimental
                    % timepoint
                        [M, I] = min(timeDifferenceAllPairs,[],1);
                    % Record the paired timepoints
                        sortedPairs = nan([length(domainSizeInUM_Exp_2Fit),2]);
                        sortedPairs(:,1) = I';
                        sortedPairs(:,2) = (1:length(domainSizeInUM_Exp_2Fit))';
        
                % Pull out only the simulation indixes matching the experimental
                % timepoints
                    timeInMin_Sim_2Compare = timeInMin_Sim(sortedPairs(:,1));
                    domainSizeInUM_Sim_2Compare = domainSizeInMicrons_All(:,sortedPairs(:,1));
                    domainSizeInUMSD_Sim_2Compare = domainSizeInMicronsSD_All(:,sortedPairs(:,1));
                    domainSizeInUMCV_Sim_2Compare = domainSizeInMicronsCV_All(:,sortedPairs(:,1));
        
                % Record the different between the simulations and the experiments
                    differenceInDomainSize_OverTime = domainSizeInUM_Sim_2Compare-domainSizeInUM_Exp_2Fit';
                    differenceInDomainSizeSD_OverTime = domainSizeInUMSD_Sim_2Compare-domainSizeInUMSD_Exp_2Fit';
                    differenceInDomainSizeCV_OverTime = domainSizeInUMCV_Sim_2Compare-domainSizeInUMCV_Exp_2Fit';
            else
                'Time vectors are not equal'
                return
            end


        % Average the simulation results across simulation replicates
            % Find the unique parameter sets
                uniqueERatios = unique(ERatio);
                numUniqueERatios = length(uniqueERatios);
                uniqueVRatios = unique(VRatio);
                numUniqueVRatios = length(uniqueVRatios);        
            % Preallocate space to store the results
                costVals_RepAvg = nan(numUniqueVRatios,numUniqueERatios);
                costValsSD_RepAvg = nan(numUniqueVRatios,numUniqueERatios);
                costValsCV_RepAvg = nan(numUniqueVRatios,numUniqueERatios);  
            % Loop through each combination and calculate the replicate-averaged
            % residuals
            for ERatioNum = 1:numUniqueERatios
                for VRatioNum = 1:numUniqueVRatios
                    % Pull out the mat files matching this combination
                        matFiles2Plot = find((ERatio==uniqueERatios(ERatioNum))&(VRatio==uniqueVRatios(VRatioNum)));
                    % Find the average residuals across replicates
                        residuals = mean(differenceInDomainSize_OverTime(matFiles2Plot,:),1,'omitnan');
                        means = mean(domainSizeInUM_Sim_2Compare(matFiles2Plot,:),1,'omitnan');
                        residualsSD = mean(differenceInDomainSizeSD_OverTime(matFiles2Plot,:),1,'omitnan');
                        meanSD = mean(domainSizeInUMSD_Sim_2Compare(matFiles2Plot,:),1,'omitnan');
                        residualsCV = mean(differenceInDomainSizeCV_OverTime(matFiles2Plot,:),1,'omitnan');
                        meanCV = mean(domainSizeInUMCV_Sim_2Compare(matFiles2Plot,:),1,'omitnan');
                    % Sum the residuals
                    if sum(isnan(residuals))>=(length(residuals)-1)
                      costVals_RepAvg(VRatioNum,ERatioNum) = nan;
                      costValsSD_RepAvg(VRatioNum,ERatioNum) = nan;
                      costValsCV_RepAvg(VRatioNum,ERatioNum) = nan;
                    else
                      costVals_RepAvg(VRatioNum,ERatioNum) = sum(abs(residuals)./means,'omitnan');
                      costValsSD_RepAvg(VRatioNum,ERatioNum) = sum(abs(residualsSD)./meanSD,'omitnan');
                      costValsCV_RepAvg(VRatioNum,ERatioNum) = sum(abs(residualsCV)./meanCV,'omitnan');
                    end
               end
            end

        % Plot the mean domain size at 48 hours as a function of the parameters 
            figure(1)
            subplot(2,4,1)
            imagesc(uniqueERatios,uniqueVRatios,log10(sortState_48hr_heatmap))
            title('Domain size at 48 hrs')
            c = colorbar;
            c.Label.String = 'Domain size (um)';
            c.TicksMode='manual';
            cTicks = log10(2.^(0:0.5:9));
            c.Ticks = cTicks;
            cTickLabels = arrayfun(@num2str, int16(10.^cTicks),'UniformOutput', 0);
            c.TickLabels = cTickLabels;
            set(gca,'YScale','log')
            xlabel(' kT_{eff}/E_{homo} (unitless)')
            ylabel(' E_{homo}/viscosity')
            set(gca,'YDir','normal') 
            set(gca,'FontName','Helvetica','FontSize',5); 
        
        
        % Plot the mean domain size at 48 hours as a function of the parameters 
            figure(1)
            subplot(2,4,2)
            imagesc(uniqueERatios,uniqueVRatios,log10(sortStateSD_48hr_heatmap))
            title('Domain size SD at 48 hrs')
            c = colorbar;
            c.Label.String = 'Domain size SD (um)';
            c.TicksMode='manual';
            cTicks = log10(2.^(0:0.5:9));
            c.Ticks = cTicks;
            cTickLabels = arrayfun(@num2str, int16(10.^cTicks),'UniformOutput', 0);
            c.TickLabels = cTickLabels;
            set(gca,'YScale','log')
            xlabel(' kT_{eff}/E_{homo} (unitless)')
            ylabel(' E_{homo}/viscosity')
            set(gca,'YDir','normal') 
            set(gca,'FontName','Helvetica','FontSize',5); 
        
        % Plot the CV of the domain size at 48 hours as a function of the parameters 
            figure(1)
            subplot(2,4,3)
            imagesc(uniqueERatios,uniqueVRatios,sortStateCV_48hr_heatmap)
            title('Domain size CV at 48 hrs')
            c = colorbar;
            c.Label.String = 'Domain size CV (um)';
            set(gca,'YScale','log')
            xlabel(' kT_{eff}/E_{homo} (unitless)')
            ylabel(' E_{homo}/viscosity')
            set(gca,'YDir','normal') 
            set(gca,'FontName','Helvetica','FontSize',5); 

        % Plot the fluidity for each set of parameters
            subplot(2,4,4)
            imagesc(uniqueERatios,uniqueVRatios,log10(fluidity_min),'alphadata',alpha)
            title('Fluidity (50% sorted)')
            c = colorbar;
            c.Label.String = 'log_{10}(Fluidity, min^{-1})';
            colormap(gca,jet*0.9)
            % c.TicksMode='manual';
            % cTicks = log10(2.^(0:0.5:9));
            % c.Ticks = cTicks;
            % cTickLabels = arrayfun(@num2str, int16(10.^cTicks),'UniformOutput', 0);
            % c.TickLabels = cTickLabels;
            set(gca,'YScale','log')
            xlabel('kT_{eff}/E_{homo} (unitless)')
            ylabel('E_{homo}/viscosity')
            set(gca,'YDir','normal') 
            set(gca,'FontName','Helvetica','FontSize',5); 

        % Plot a heatmap of the difference between the simulation and the
        % experiment across all simulation parameters
            % Mean domain size    
                figure(1)
                subplot(2,4,5)
                imagesc(uniqueERatios,uniqueVRatios,log10(costVals_RepAvg))
                title('Difference from experiment')
                c = colorbar;
                c.Label.String = 'Summed residuals (um)';
                c.TicksMode='manual';
                cTicks = log10(10.^(1:6));
                c.Ticks = cTicks;
                cTickLabels = arrayfun(@num2str, 10.^cTicks, 'UniformOutput', 0);
                c.TickLabels = cTickLabels;
                set(gca,'YScale','log')
                xlabel(' kT_{eff}/E_{homo} (unitless)')
                ylabel(' E_{homo}/viscosity (?)')
                set(gca,'YDir','normal')    
                set(gca,'FontName','Helvetica','FontSize',5);      
            % SD of domain size
                figure(1)
                subplot(2,4,6)
                imagesc(uniqueERatios,uniqueVRatios,log10(costValsSD_RepAvg))
                title('Difference from experiment')
                c = colorbar;
                c.Label.String = 'Summed residuals of SD (um)';
                c.TicksMode='manual';
                cTicks = log10(10.^(1:6));
                c.Ticks = cTicks;
                cTickLabels = arrayfun(@num2str, 10.^cTicks, 'UniformOutput', 0);
                c.TickLabels = cTickLabels;
                set(gca,'YScale','log')
                xlabel(' kT_{eff}/E_{homo} (unitless)')
                ylabel(' E_{homo}/viscosity (?)')
                set(gca,'YDir','normal') 
                set(gca,'FontName','Helvetica','FontSize',5);
            % CV of domain size
                figure(1)
                subplot(2,4,7)
                imagesc(uniqueERatios,uniqueVRatios,costValsCV_RepAvg)
                title('Difference from experiment')
                c = colorbar;
                c.Label.String = 'Summed residuals of CV (um)';
                % c.TicksMode='manual';
                % cTicks = log10(10.^(1:6));
                % c.Ticks = cTicks;
                % cTickLabels = arrayfun(@num2str, 10.^cTicks, 'UniformOutput', 0);
                % c.TickLabels = cTickLabels;
                set(gca,'YScale','log')
                xlabel(' kT_{eff}/E_{homo} (unitless)')
                ylabel(' E_{homo}/viscosity (?)')
                set(gca,'YDir','normal') 
                set(gca,'FontName','Helvetica','FontSize',5);

        % Find the best fit simulation to the experiments
    
            % Find the simuations with minimam difference
                [minVal, minIdx] = min(costVals_RepAvg(:)+1*costValsSD_RepAvg(:));
            % Convert to subscript indices
                [row col] = ind2sub(size(costVals_RepAvg),minIdx);
            % Pull out the optimum
                optERatio = uniqueERatios(col);
                optVRatio = uniqueVRatios(row);
            % Plot the optimum fit on top of the heatmap
                    figure(1)
                    subplot(2,4,1)
                    hold on;
                    scatter(uniqueERatios(col),uniqueVRatios(row),'w*')
                    set(gca,'FontName','Helvetica','FontSize',5); 
                    hold off;
                    subplot(2,4,2)
                    hold on;
                    scatter(uniqueERatios(col),uniqueVRatios(row),'w*')
                    set(gca,'FontName','Helvetica','FontSize',5); 
                    hold off;
                    subplot(2,4,3)
                    hold on;
                    scatter(uniqueERatios(col),uniqueVRatios(row),'w*')
                    set(gca,'FontName','Helvetica','FontSize',5); 
                    hold off;
                    subplot(2,4,4)
                    hold on;
                    scatter(uniqueERatios(col),uniqueVRatios(row),'w*')
                    set(gca,'FontName','Helvetica','FontSize',5); 
                    hold off;
                    subplot(2,4,5)
                    hold on;
                    scatter(uniqueERatios(col),uniqueVRatios(row),'w*')
                    set(gca,'FontName','Helvetica','FontSize',5); 
                    hold off;
                    subplot(2,4,6)
                    hold on;
                    scatter(uniqueERatios(col),uniqueVRatios(row),'w*')
                    set(gca,'FontName','Helvetica','FontSize',5); 
                    hold off;
                    subplot(2,4,7)
                    hold on;
                    scatter(uniqueERatios(col),uniqueVRatios(row),'w*')
                    hold off;
                    set(gca,'FontName','Helvetica','FontSize',5); 


        % Save a figure of the best fit
            % Choose the video file path
                figFilePath = [figOutputFolderPath 'Fig6_PhaseSpace_' dataNames{iSeries}];
            % Choose the figure size and position
                figHandle = figure(1);
                figHandle.Position =  [145   55   975   250];   
            % Save the file as a pdf
                exportgraphics(gcf,[figFilePath '.pdf'],'ContentType','vector')
            % Save as a png                
                saveas(gcf,[figFilePath '.png'])

        % Plot a comparison of the sorting behavior

            % Initialize the legend text
                legendText={};
            
            % Initialize the plot counter
                plotNum=2;
            
            if sum(eqVec)==length(eqVec)
            
                % Subset the experimental data to times available in the simulation
                    timeInMin_Exp_2Fit = timeInMin_Exp(timeInMin_Exp<=max(timeInMin_Sim));
                    domainSizeInUM_Exp_2Fit = domainSizeInUM_Exp(timeInMin_Exp<=max(timeInMin_Sim));
                    domainSizeInUMSD_Exp_2Fit = domainSizeInUMSD_Exp(timeInMin_Exp<=max(timeInMin_Sim));
            
                % Find all replicates that match the optimum
                    matFiles2Plot = find((optERatio==ERatio)&...
                        (optVRatio==VRatio));

                minDomainSize = Inf;
                maxDomainSize = 0;
                minDomainSizeSD = Inf;
                maxDomainSizeSD = 0;
            
                % Loop through each replicate and plot the comparison    
                for repNum=1:length(matFiles2Plot)
            
                    % Pull out the simulation information
                        domainSizeInUM_Sim_ThisRep = domainSizeInMicrons_All(matFiles2Plot(repNum),:);
                        domainSizeInUMSD_Sim_ThisRep = domainSizeInMicronsSD_All(matFiles2Plot(repNum),:);
                        domainSizeInUM_Sim_2Compare_ThisRep = domainSizeInUM_Sim_2Compare(matFiles2Plot(repNum),:)';
                        domainSizeInUMSD_Sim_2Compare_ThisRep = domainSizeInUMSD_Sim_2Compare(matFiles2Plot(repNum),:)';
            
                        timeInMin_Sim = timeInMin_All(matFileNum,:);
                
                    % Plot the simulation at full resolution
                        figure(2)
                        subplot(2,5,1)
                        if repNum==1
                            plot(timeInMin_Sim./60,domainSizeInUM_Sim_ThisRep,'k-')
                            hold on;
                        else
                            plot(timeInMin_Sim./60,domainSizeInUM_Sim_ThisRep,'k-','HandleVisibility','off')
                            hold on;
                        end
                        subplot(2,5,6)
                        if repNum==1
                            plot(timeInMin_Sim./60,domainSizeInUMSD_Sim_ThisRep,'k-')
                            hold on;
                        else
                            plot(timeInMin_Sim./60,domainSizeInUMSD_Sim_ThisRep,'k-','HandleVisibility','off')
                            hold on;
                        end
                

                    % Record the minimum so far
                        minDomainSize = min([minDomainSize,min(domainSizeInUM_Sim_ThisRep)]);
                        maxDomainSize = max([maxDomainSize,max(domainSizeInUM_Sim_ThisRep)]);
                        minDomainSizeSD = min([minDomainSizeSD,min(domainSizeInUMSD_Sim_ThisRep)]);
                        maxDomainSizeSD = max([maxDomainSizeSD,max(domainSizeInUMSD_Sim_ThisRep)]);

                    % Update the plot counter
                        plotNum = plotNum + 1;
            
                end
            
            end
            
            % Record the minimum so far
                minDomainSize = min([minDomainSize,min(domainSizeInUM_Exp)]);
                maxDomainSize = max([maxDomainSize,max(domainSizeInUM_Exp)]);
                minDomainSizeSD = min([minDomainSizeSD,min(domainSizeInUMSD_Exp)]);
                maxDomainSizeSD = max([maxDomainSizeSD,max(domainSizeInUMSD_Exp)]);

            % Plot the experimental data
                subplot(2,5,1)
                scatter(timeInMin_Exp./60,domainSizeInUM_Exp,3,'ro','filled')
                % Clean up the plot
                    hold off;
                    xlim([0 max(timeInMin_Exp)]./60)
                    ylim([0 85])
                    ylim([minDomainSize maxDomainSize])
                    xlabel('Time (hrs)')
                    ylabel('Mean domain size (um)')
                    % Print the best fit parameters to the legend
                        legendText{1} = sprintf('kT = %1.2f*E_{homo} \n E_{homo} = %1.0e*viscosity \n %i replicates',...
                            ERatio(matFiles2Plot(repNum)), VRatio(matFiles2Plot(repNum)),length(matFiles2Plot));
                        legendText{2} = 'Experimental data';
                           % legend(legendText,'Location','Northwest')
                    set(gca,'FontName','Helvetica','FontSize',5); 
            
                subplot(2,5,6)
                scatter(timeInMin_Exp./60,domainSizeInUMSD_Exp,3,'ro','filled')    
                % Clean up the plot
                    hold off;
                    xlim([0 max(timeInMin_Exp)]./60)
                    ylim([0 85])
                    ylim([minDomainSizeSD maxDomainSizeSD])
                    xlabel('Time (hrs)')
                    ylabel('SD domain size (um)')
                    % Print the best fit parameters to the legend
                        legendText{1} = sprintf('kT = %1.2f*E_{homo} \n E_{homo} = %1.0e*viscosity \n %i replicates',...
                            ERatio(matFiles2Plot(repNum)), VRatio(matFiles2Plot(repNum)),length(matFiles2Plot));
                        legendText{2} = 'Experimental data';
                            %legend(legendText,'Location','Northwest')
                    set(gca,'FontName','Helvetica','FontSize',5); 
            
                subplot(2,5,1)
                hold off;
                xlabel('Time (hrs)')
                ylabel('Domain size (um)')
                subplot(2,5,6)
                hold off;
                xlabel('Time (hrs)')
                ylabel('Domain size SD (um)')

        % Save the information
            ERatio_BestFitExp(seriesCounter) = optERatio;
            VRatio_BestFitExp(seriesCounter) = optVRatio;
            domainSizeInMicrons_BestFitExp{seriesCounter} = domainSizeInUM_Exp;
            domainSizeInMicronsSD_BestFitExp{seriesCounter} = domainSizeInUMSD_Exp;
            timeInMin_BestFitExp{seriesCounter} = timeInMin_Exp;


        % Plot the comparison of the time lapse montage

            % Load the experimental data
                % Load the cell type data
                    clear cellType_OverTime
                    load(expMatFilePath,'cellType_OverTime');        
                % Pull out the file path for experimental image data
                    prefixUS = find(expMatFilePath=='_',2,'last');
                    prefixUS = prefixUS(1);
                    expImagePath =[expMatFilePath(1:prefixUS) 'MaxIP_MATLAB.tiff'];

            % Pull out the simulated data for the best fit 
                % Pull out the file path
                    matFilePath = [simDataFolderPath matFiles(matFiles2Plot(repNum)).name];
                % Load the relevant data
                    load(matFilePath,'cellID_onGrid_OverTime','cellType_onGrid_OverTime','timeInMin_OverTime')

            % Pull out the timepoints to plot in the montage

                % Experiments

                    % Choose the time points to plot in min
                        timeInMin2Plot_Ideal = [0 30 3*60 (60*18)]; 
        
                    % Find the best matching timepoints in the experimental data
        
                    % Find the associated time in the simulations that are closest to the
                    % timepoints selected for plotting
                        % Calculate the time different between all possible pairs of
                        % the experimental time and the actual timepoints
                            timeDifferenceAllPairs_Exp = abs(timeInMin_Exp_2Fit(:)-timeInMin2Plot_Ideal(:).');
                        % Find the closest simulated timepoint to each experimental
                        % timepoint
                            [M, timePointIndex_Exp] = min(timeDifferenceAllPairs_Exp,[],1);
                        % Pull out the actual time values
                            timeInMin2Plot_Exp = timeInMin_Exp_2Fit(timePointIndex_Exp);

                % Simulations

                    % Find the associated time in the simulations that are closest to the
                    % timepoints selected for plotting
                        % Calculate the time different between all possible pairs of
                        % the experimental time and the actual timepoints
                            timeDifferenceAllPairs_Sim = abs(timeInMin_OverTime(:)-timeInMin2Plot_Ideal(:).');
                        % Find the closest simulated timepoint to each experimental
                        % timepoint
                            [M, timePointIndex_Sim] = min(timeDifferenceAllPairs_Sim,[],1);
                        % Pull out the actual time values
                            timeInMin2Plot_Sim = timeInMin_OverTime(timePointIndex_Sim);

            % Choose the x limits
                xlimMax = max([size(cellType_OverTime,1)*micronsPerPixel_X,...
                    size(cellID_onGrid_OverTime,1)*cellRadius*2/1000]);

            % Read the file (but do not load the images)
                reader_MIP = bfGetReader(expImagePath);

           % Initialize the counter
                tpCounter = 0;

                figure(2)
            for numTP = timePointIndex_Exp

                % Update the counter
                    tpCounter = tpCounter+1;

                % Load the image for the experiments
                    % Pull out the appropriate planes
                        iPlane1 = reader_MIP.getIndex(1 - 1, 1 -1, numTP - 1) + 1;
                        iPlane2 = reader_MIP.getIndex(1 - 1, 2 -1, numTP - 1) + 1;
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

                % Create an associated RGB image
                    MIP_RGB = nan(size(cellType1ChannelIdxImage_MIP,1),size(cellType1ChannelIdxImage_MIP,2),3);
                    MIP_RGB(:,:,2) = imadjust(cellType1ChannelIdxImage_MIP);
                    MIP_RGB(:,:,1) = imadjust(cellType2ChannelIdxImage_MIP);
                    MIP_RGB(:,:,3) = imadjust(cellType2ChannelIdxImage_MIP);


                % Plot the experiments
                    % Plot the MIP
                        subplot(2,5,1+tpCounter)
                        imagesc(1:size(cellType_OverTime,1)*micronsPerPixel_X,...
                            1:size(cellType_OverTime,1)*micronsPerPixel_X,...
                            MIP_RGB)      
                        xlabel('Microns')
                        ylabel('Microns')
                        xlim([0 xlimMax])
                        ylim([0 xlimMax])
                        title('Experiment MIP')
                        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');                        
                        if timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter))==0
                        title(sprintf('Time = %1.1f min',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))) 
                        elseif abs(timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))<0.4
                        title(sprintf('Time = %1.1f sec',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter))*60)) 
                        elseif abs(timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))<59
                        title(sprintf('Time = %1.1f min',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))) 
                        else
                        title(sprintf('Time = %1.1f hrs',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter))/60)) 
                        end
                        set(gca,'FontName','Helvetica','FontSize',5); 
                        
                    % Plot the cell type
                        subplot(2,5,5+1+tpCounter)
                        imagesc(1:size(cellID_onGrid_OverTime,1)*cellRadius*2/1000,...
                            1:size(cellID_onGrid_OverTime,1)*cellRadius*2/1000,...
                            squeeze(cellType_onGrid_OverTime(:,:,timePointIndex_Sim(tpCounter))))                        
                        if timeInMin_OverTime(timePointIndex_Sim(tpCounter))==0
                        title(sprintf('Time = %1.1f min',timeInMin_OverTime(timePointIndex_Sim(tpCounter)))) 
                        elseif abs(timeInMin_OverTime(timePointIndex_Sim(tpCounter)))<0.4
                        title(sprintf('Time = %1.1f sec',timeInMin_OverTime(timePointIndex_Sim(tpCounter))*60)) 
                        elseif abs(timeInMin_OverTime(timePointIndex_Sim(tpCounter)))<59
                        title(sprintf('Time = %1.1f min',timeInMin_OverTime(timePointIndex_Sim(tpCounter)))) 
                        else
                        title(sprintf('Time = %1.1f hrs',timeInMin_OverTime(timePointIndex_Sim(tpCounter))/60)) 
                        end
                       % set(gca,'Color','k')
                        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
                        xlim([0 xlimMax])
                        ylim([0 xlimMax])
                        text(10,1000,legendText{1},'FontName','Helvetica','FontSize',5); 
                        drawnow
                        set(gca,'FontName','Helvetica','FontSize',5); 


            end



        % Save a figure of the best fit
            % Pull out the file path
                figFilePath = [figOutputFolderPath 'Fig6_TLM_' ...
                    dataNames{iSeries} '_' matFiles(matFiles2Plot(repNum)).name];
            % Choose the video file path
                figFilePath = [figFilePath(1:(find(figFilePath=='.',1,'last'))-1)];
            % Choose the figure size and position
                figHandle = figure(2);
                figHandle.Position =  [175  400   716   232];
            % Save the file as a pdf
                exportgraphics(gcf,[figFilePath '.pdf'],'ContentType','vector')
            % Save as a png                
                saveas(gcf,[figFilePath '.png'])


        % Pull out the timepoints to plot in the video

                % Experiments

                    % Choose the time points to plot in min
                        timeInMin2Plot_Ideal = timeInMin_Exp_2Fit; 
        
                    % Find the best matching timepoints in the experimental data
        
                    % Find the associated time in the simulations that are closest to the
                    % timepoints selected for plotting
                        % Calculate the time different between all possible pairs of
                        % the experimental time and the actual timepoints
                            timeDifferenceAllPairs_Exp = abs(timeInMin_Exp_2Fit(:)-timeInMin2Plot_Ideal(:).');
                        % Find the closest simulated timepoint to each experimental
                        % timepoint
                            [M, timePointIndex_Exp] = min(timeDifferenceAllPairs_Exp,[],1);
                        % Pull out the actual time values
                            timeInMin2Plot_Exp = timeInMin_Exp_2Fit(timePointIndex_Exp);

                % Simulations

                    % Find the associated time in the simulations that are closest to the
                    % timepoints selected for plotting
                        % Calculate the time different between all possible pairs of
                        % the experimental time and the actual timepoints
                            timeDifferenceAllPairs_Sim = abs(timeInMin_OverTime(:)-timeInMin2Plot_Ideal(:).');
                        % Find the closest simulated timepoint to each experimental
                        % timepoint
                            [M, timePointIndex_Sim] = min(timeDifferenceAllPairs_Sim,[],1);
                        % Pull out the actual time values
                            timeInMin2Plot_Sim = timeInMin_OverTime(timePointIndex_Sim);

            % Choose the x limits
                xlimMax = max([size(cellType_OverTime,1)*micronsPerPixel_X,...
                    size(cellID_onGrid_OverTime,1)*cellRadius*2/1000]);

            % Read the file (but do not load the images)
                reader_MIP = bfGetReader(expImagePath);

           % Initialize the counter
                tpCounter = 0;

            % Choose and create a directory to save the movie

                % Pull out the file path
                    videoFilePath = [figOutputFolderPath matFiles(matFiles2Plot(repNum)).name];
                % Choose the video file path
                    videoFilePath = [videoFilePath(1:(find(videoFilePath=='.',1,'last'))-1) '_' dataNames{iSeries} '.avi'];
                % Start the video writer
                     v = VideoWriter(videoFilePath,'Uncompressed AVI');
                     v.FrameRate = round(7);
                     open(v);

                figure(30)
            for numTP = timePointIndex_Exp

                % Update the counter
                    tpCounter = tpCounter+1;

                % Load the image for the experiments
                    % Pull out the appropriate planes
                        iPlane1 = reader_MIP.getIndex(1 - 1, 1 -1, numTP - 1) + 1;
                        iPlane2 = reader_MIP.getIndex(1 - 1, 2 -1, numTP - 1) + 1;
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

                % Create an associated RGB image
                    MIP_RGB = nan(size(cellType1ChannelIdxImage_MIP,1),size(cellType1ChannelIdxImage_MIP,2),3);
                    MIP_RGB(:,:,2) = imadjust(cellType1ChannelIdxImage_MIP);
                    MIP_RGB(:,:,1) = imadjust(cellType2ChannelIdxImage_MIP);
                    MIP_RGB(:,:,3) = imadjust(cellType2ChannelIdxImage_MIP);


                % Plot the comparison
                    % Plot the experiment MIP
                        figure(30)
                        subplot(2,1,1)
                        imagesc(1:size(cellType_OverTime,1)*micronsPerPixel_X,...
                            1:size(cellType_OverTime,1)*micronsPerPixel_X,...
                            MIP_RGB)      
                        xlabel('Microns')
                        ylabel('Microns')
                        axis equal
                        xlim([0 xlimMax])
                        ylim([0 xlimMax])
                        title('Experiment MIP')
                        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');                        
                        if timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter))==0
                        title(sprintf('Time = %1.1f min',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))) 
                        elseif abs(timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))<0.4
                        title(sprintf('Time = %1.1f sec',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter))*60)) 
                        elseif abs(timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))<59
                        title(sprintf('Time = %1.1f min',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter)))) 
                        else
                        title(sprintf('Time = %1.1f hrs',timeInMin_Exp_2Fit(timePointIndex_Exp(tpCounter))/60)) 
                        end
                        set(gca,'FontName','Helvetica','FontSize',5); 
                    % Plot the simualtion
                        subplot(2,1,2)
                        imagesc(1:size(cellID_onGrid_OverTime,1)*cellRadius*2/1000,...
                            1:size(cellID_onGrid_OverTime,1)*cellRadius*2/1000,...
                            squeeze(cellType_onGrid_OverTime(:,:,timePointIndex_Sim(tpCounter))))                        
                        if timeInMin_OverTime(timePointIndex_Sim(tpCounter))==0
                        title(sprintf('Time = %1.1f min',timeInMin_OverTime(timePointIndex_Sim(tpCounter)))) 
                        elseif abs(timeInMin_OverTime(timePointIndex_Sim(tpCounter)))<0.4
                        title(sprintf('Time = %1.1f sec',timeInMin_OverTime(timePointIndex_Sim(tpCounter))*60)) 
                        elseif abs(timeInMin_OverTime(timePointIndex_Sim(tpCounter)))<59
                        title(sprintf('Time = %1.1f min',timeInMin_OverTime(timePointIndex_Sim(tpCounter)))) 
                        else
                        title(sprintf('Time = %1.1f hrs',timeInMin_OverTime(timePointIndex_Sim(tpCounter))/60)) 
                        end
                       % set(gca,'Color','k')
                        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
                        axis equal
                        xlim([0 xlimMax])
                        ylim([0 xlimMax])
                        text(10,1000,legendText{1},'FontName','Helvetica','FontSize',5); 
                        drawnow
                        set(gca,'FontName','Helvetica','FontSize',5); 

                % And save the video
                    % Choose the figure size and position
                        figHandle = figure(30);
                        figHandle.Position =  [500   400   270   350];
                    % Grab the figure frame
                        F = getframe(figHandle);
                    % Append the fram to the video
                        writeVideo(v,F);

            end

                close(v);

    end







end

%% Fig. 6a-b -- Plot a figure of the results

cmap = [1 0 0; 0 1 1];
markers = {'o','s','d','^'};
filled = {'filled',''};

domainSizeInMicrons_BestFitExp_All = nan([length(Day_BestFitExp),length(timeInMin_BestFitExp{1})]);
domainSizeInMicronsSD_BestFitExp_All = nan([length(Day_BestFitExp),length(timeInMin_BestFitExp{1})]);
domainSizeInMicronsCV_BestFitExp_All = nan([length(Day_BestFitExp),length(timeInMin_BestFitExp{1})]);

for seriesNum = 1:length(Day_BestFitExp)

    % Plot the experimental data
        figure(3)
        subplot(1,3,1)
        scatter(timeInMin_BestFitExp{seriesNum}./60,domainSizeInMicrons_BestFitExp{seriesNum},...
            3,markers{ReplicateWithinDay_BestFitExp(seriesNum)+2*(Day_BestFitExp(seriesNum)-1)},'filled','MarkerFaceColor',...
            cmap(CadCombo_BestFitExp(seriesNum),:)./ExpressionLevel_BestFitExp(seriesNum), ...
            'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
        hold on;
        subplot(1,3,2)
        scatter(timeInMin_BestFitExp{seriesNum}./60,domainSizeInMicronsSD_BestFitExp{seriesNum},...
            3,markers{ReplicateWithinDay_BestFitExp(seriesNum)+2*(Day_BestFitExp(seriesNum)-1)},'filled','MarkerFaceColor',...
            cmap(CadCombo_BestFitExp(seriesNum),:)./ExpressionLevel_BestFitExp(seriesNum), ...
            'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
        hold on;
        subplot(1,3,3)
        scatter(timeInMin_BestFitExp{seriesNum}./60,domainSizeInMicronsSD_BestFitExp{seriesNum}./domainSizeInMicrons_BestFitExp{seriesNum},...
            3,markers{ReplicateWithinDay_BestFitExp(seriesNum)+2*(Day_BestFitExp(seriesNum)-1)},'filled','MarkerFaceColor',...
            cmap(CadCombo_BestFitExp(seriesNum),:)./ExpressionLevel_BestFitExp(seriesNum), ...
            'MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
        hold on;
            
        % Save the results
            domainSizeInMicrons_BestFitExp_All(seriesNum,:) = domainSizeInMicrons_BestFitExp{seriesNum};
            domainSizeInMicronsSD_BestFitExp_All(seriesNum,:) = domainSizeInMicronsSD_BestFitExp{seriesNum};
            domainSizeInMicronsCV_BestFitExp_All(seriesNum,:) = domainSizeInMicronsSD_BestFitExp{seriesNum}./domainSizeInMicrons_BestFitExp{seriesNum};

end

    % Label the plot
        subplot(1,3,1)
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==1),:)),'-','color',...
           cmap(1,:)./1,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==1),:)),...
            3,'o','filled','MarkerFaceColor',cmap(1,:)./1,'MarkerEdgeColor',[0 0 0])        
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==2),:)),'-','color',...
           cmap(1,:)./2,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==2),:)),...
            3,'o','filled','MarkerFaceColor',cmap(1,:)./2,'MarkerEdgeColor',[0 0 0])        
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==1),:)),'-','color',...
           cmap(2,:)./1,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==1),:)),...
            3,'o','filled','MarkerFaceColor',cmap(2,:)./1,'MarkerEdgeColor',[0 0 0])      
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==2),:)),'-','color',...
           cmap(2,:)./2,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicrons_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==2),:)),...
            3,'o','filled','MarkerFaceColor',cmap(2,:)./2,'MarkerEdgeColor',[0 0 0])
        xlim([0 18])
        hold off;
        xlabel('Time (hrs)')
        ylabel('Domain size (um)')
        set(gca,'FontName','Helvetica','FontSize',5); 
        box on;
        subplot(1,3,2)
        xlim([0 18])
        hold off;
        xlabel('Time (hrs)')
        ylabel('Domain size SD (um)')
        set(gca,'FontName','Helvetica','FontSize',5); 
        box on;
        subplot(1,3,3)
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==1),:)),'-','color',...
           cmap(1,:)./1,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==1),:)),...
            3,'o','filled','MarkerFaceColor',cmap(1,:)./1,'MarkerEdgeColor',[0 0 0])        
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==2),:)),'-','color',...
           cmap(1,:)./2,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==2),:)),...
            3,'o','filled','MarkerFaceColor',cmap(1,:)./2,'MarkerEdgeColor',[0 0 0])        
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==1),:)),'-','color',...
           cmap(2,:)./1,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==1),:)),...
            3,'o','filled','MarkerFaceColor',cmap(2,:)./1,'MarkerEdgeColor',[0 0 0])      
        plot(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==2),:)),'-','color',...
           cmap(2,:)./2,'LineWidth',1)
        scatter(timeInMin_BestFitExp{1}./60,mean(domainSizeInMicronsCV_BestFitExp_All((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==2),:)),...
            3,'o','filled','MarkerFaceColor',cmap(2,:)./2,'MarkerEdgeColor',[0 0 0])
        xlim([0 18])
        ylim([0.3 0.4])
        hold off;
        xlabel('Time (hrs)')
        ylabel('Domain size CV (unitless)')
        set(gca,'FontName','Helvetica','FontSize',5); 
        box on;

        % Save a figure of the best fit
            % Pull out the file path
                figFilePath = [figOutputFolderPath 'Fig6ab_ExpmtSortingByCondition'];
            % Choose the figure size and position
                figHandle = figure(3);
                figHandle.Position =  [330  575   315   90];
            % Save the file as a pdf
                exportgraphics(gcf,[figFilePath '.pdf'],'ContentType','vector')
            % Save as a png                
                saveas(gcf,[figFilePath '.png'])



%% Fig. 6d -- Plot the best fit parameters for each condition on top of the heatmap 

tIdx_18hrs = find(timeInMin_All(1,:)>(18*60),1,'first');

% Calcluate the mean, SD, and CV of the domain size at 18 hrs for all
% simulation parameters
    domainSizeInMicrons_18hr = nanmean(domainSizeInMicrons_All(:,(tIdx_18hrs-50):tIdx_18hrs),2);
    domainSizeInMicronsSD_18hr = nanmean(domainSizeInMicronsSD_All(:,(tIdx_18hrs-50):tIdx_18hrs),2);
    domainSizeInMicronsCV_18hr = nanmean(domainSizeInMicronsCV_All(:,(tIdx_18hrs-50):tIdx_18hrs),2);

% Create a matrix to store the domain sizes
    sortState_18hr_heatmap = nan([numUniqueVRatios,numUniqueERatios]);
    sortStateSD_18hr_heatmap = nan([numUniqueVRatios,numUniqueERatios]);
    sortStateCV_18hr_heatmap = nan([numUniqueVRatios,numUniqueERatios]);

% Loop through each set of parameters and record the domain size
for ERatioNum = 1:numUniqueERatios
    for VRatioNum = 1:numUniqueVRatios
       
        % Pull out the mat files matching this combination
            matFiles2Plot = find((ERatio==uniqueERatios(ERatioNum))&(VRatio==uniqueVRatios(VRatioNum)));

        % Fill in the heatmap
        if any(~isnan(domainSizeInMicrons_18hr(matFiles2Plot)))
            sortState_18hr_heatmap(VRatioNum,ERatioNum) = mean(domainSizeInMicrons_18hr(matFiles2Plot),'omitnan');
            sortStateSD_18hr_heatmap(VRatioNum,ERatioNum) = mean(domainSizeInMicronsSD_18hr(matFiles2Plot),'omitnan');
            sortStateCV_18hr_heatmap(VRatioNum,ERatioNum) = mean(domainSizeInMicronsCV_18hr(matFiles2Plot),'omitnan');
        end
    end
end


% Plot the mean domain size at 48 hours as a function of the parameters 
    figure(4)
    imagesc(uniqueERatios,uniqueVRatios*10,log10(sortState_18hr_heatmap))
    hold on;
    title('Domain size at 18 hrs')
    c = colorbar;
    %colormap(gca,slanCM('gem')*0.9)
    colormap(gca,gray)
        c.Label.String = 'Domain size (um)';
    c.TicksMode='manual';
    cTicks = log10(2.^(0:0.5:9));
    c.Ticks = cTicks;
    cTickLabels = arrayfun(@num2str, int16(10.^cTicks),'UniformOutput', 0);
    c.TickLabels = cTickLabels;
    %clim([1 2.7])
    set(gca,'YScale','log')
    xlabel(' E_M (multiples of E_{homo})')
    ylabel(' E_{homo}(assuming n = 10 Pa s)')
    set(gca,'YDir','normal') 
    set(gca,'FontName','Helvetica','FontSize',5); 


cmap = [1 0 0; 0 1 1];
markers = {'o','s','d','^'};

domainSizeInMicrons_BestFitExp_All = nan([length(Day_BestFitExp),length(timeInMin_BestFitExp{1})]);
domainSizeInMicronsSD_BestFitExp_All = nan([length(Day_BestFitExp),length(timeInMin_BestFitExp{1})]);
domainSizeInMicronsCV_BestFitExp_All = nan([length(Day_BestFitExp),length(timeInMin_BestFitExp{1})]);

for seriesNum = 1:length(Day_BestFitExp)

    % Plot the experimental data
        figure(4)
        scatter(ERatio_BestFitExp(seriesNum),VRatio_BestFitExp(seriesNum)*10,...
            6,markers{ReplicateWithinDay_BestFitExp(seriesNum)+2*(Day_BestFitExp(seriesNum)-1)},'filled','MarkerFaceColor',...
            cmap(CadCombo_BestFitExp(seriesNum),:)./ExpressionLevel_BestFitExp(seriesNum), ...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        hold on;

end

    % Label the plot
        scatter(mean(ERatio_BestFitExp((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==1),:)),10*mean(VRatio_BestFitExp((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==1),:)),...
            6,'o','filled','MarkerFaceColor',cmap(1,:)./1,'MarkerEdgeColor',[1 1 1])        
        scatter(mean(ERatio_BestFitExp((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==2),:)),10*mean(VRatio_BestFitExp((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==2),:)),...
            6,'o','filled','MarkerFaceColor',cmap(1,:)./2,'MarkerEdgeColor',[1 1 1])        
        scatter(mean(ERatio_BestFitExp((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==1),:)),10*mean(VRatio_BestFitExp((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==1),:)),...
            6,'o','filled','MarkerFaceColor',cmap(2,:)./1,'MarkerEdgeColor',[1 1 1])      
        scatter(mean(ERatio_BestFitExp((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==2),:)),10*mean(VRatio_BestFitExp((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==2),:)),...
            6,'o','filled','MarkerFaceColor',cmap(2,:)./2,'MarkerEdgeColor',[1 1 1])
        hold off;
        set(gca,'ytick',[10^3 10^4 10^5 10^6 10^7 10^8])
        set(gca,'xtick',[0 0.25 0.5 0.75 1 1.25 1.5])
    
    % Save a figure of the best fit
        % Pull out the file path
            figFilePath = [figOutputFolderPath 'Fig6d_ExpmtFitOnHeatmap'];
        % Choose the figure size and position
            figHandle = figure(4);
            figHandle.Position =  [488  550   160   110];
        % Save the file as a pdf
            exportgraphics(gcf,[figFilePath '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[figFilePath '.png'])


%% Fig. 6e -- Plot and fit how motility scales with the adhesion energy

for seriesNum = 1:length(Day_BestFitExp)

    % Plot the experimental data
        figure(5)
        scatter(10*VRatio_BestFitExp(seriesNum),ERatio_BestFitExp(seriesNum).*10*VRatio_BestFitExp(seriesNum),...
            60,markers{ReplicateWithinDay_BestFitExp(seriesNum)+2*(Day_BestFitExp(seriesNum)-1)},'filled','MarkerFaceColor',...
            cmap(CadCombo_BestFitExp(seriesNum),:)./ExpressionLevel_BestFitExp(seriesNum), ...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        hold on;

end

set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'xtick',[10^3 10^4 10^5 10^6 10^7 10^8])
set(gca,'ytick',[10^3 10^4 10^5 10^6 10^7 10^8])
xlabel(' E_{M} (multiples of k_BT)')
ylabel(' E_{homo} (multiples of k_BT)')
xlim([2*10^4 2*10^6])
ylim([2*10^4 2*10^6])

% Prepare the data to fit
    %newY = ERatio_BestFitExp.*VRatio_BestFitExp;
    %newLogY = log(newY);
% Set up fittype and options.
    ft2 = fittype('(m*x)+b');
% Fit model to data.
    [fitresult2, gof] = fit( log(10*VRatio_BestFitExp),log(10*VRatio_BestFitExp.*ERatio_BestFitExp),   ft2 );
% Plot fit with data.
    plot(sort(10*VRatio_BestFitExp),exp(fitresult2.b)*(sort(10*VRatio_BestFitExp).^fitresult2.m),'k--');


        x2 = [10.^([1.5:8.5]), fliplr(10.^([1.5:8.5]))];
        inBetween = [10.^([1.5:8.5])*5, fliplr(10.^([1.5:8.5])/1.5)];
        fill(x2, inBetween, 'k','FaceAlpha',0.1,'EdgeAlpha',0);

    hold off;
% Label the plot
    legend('Expmt best fit',...
        sprintf('kT ∝ E_{homo}^{%1.2f}',fitresult2.m),'Narrow sorting window',...
        'Location','Northwest')
    box on;
    % ylabel(' kT_{eff}/viscosity')
    % xlabel(' E_{homo}/viscosity')



    % Save a figure of the best fit
        % Pull out the file path
            figFilePath = [figOutputFolderPath 'Fig6e_ExpmtFitScalekTwithE'];
        % Choose the figure size and position
            figHandle = figure(5);
            figHandle.Position =  [488  438   250   220];
        % Save the file as a pdf
            exportgraphics(gcf,[figFilePath '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[figFilePath '.png'])


%% Calculate the mean intensity for each experimental condition

% Create a space to store the best fit values
    cellType_All = cell(16,1);
    rawImg_CellType1_All = cell(16,1);
    rawImg_CellType2_All = cell(16,1);

% Create a space to store the best fit values
    CadCombo_BestFitExp = nan(16,1);
    Day_BestFitExp = nan(16,1);
    ReplicateWithinDay_BestFitExp = nan(16,1);
    ExpressionLevel_BestFitExp = nan(16,1);

% Create a counter
    seriesCounter = 0;

% Loop through each dataset
for datasetNum = 1:2

    % Select the data to fit
    if datasetNum==1
        % Input the image file path
            imageFilePath = [expmtMasterFolderPath 'Processed_Data\240418_spheroid_assay\L929_2h_220um_timeseries.nd2'];
        % Input the time delay between plating and imaging
            timeDelayInMin = -120;
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
        % Choose the dataset names
            dataNames = strcat(repmat({dateVal},length(seriesNums),1),...
                repmat({'_S'},length(seriesNums),1),...
                arrayfun(@num2str, seriesNums', 'UniformOutput', 0),...
                repmat({'_'},length(seriesNums),1),...
                arrayfun(@num2str, numCells(seriesNums)', 'UniformOutput', 0)',...
                repmat({'_'},length(seriesNums),1),...
                cellType1Name(seriesNums)',...
                repmat({'_'},length(seriesNums),1),...
                cellType2Name(seriesNums)');
    elseif datasetNum==2
        % Input the image file path
            imageFilePath = [expmtMasterFolderPath 'Processed_Data\240420_spheroid_assay\L929_2h_220um_timeseries.nd2'];
        % Input the time delay between plating and imaging
            timeDelayInMin = -120;
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
        % Choose the dataset names
            dataNames = strcat(repmat({dateVal},length(seriesNums),1),...
                repmat({'_S'},length(seriesNums),1),...
                arrayfun(@num2str, seriesNums', 'UniformOutput', 0),...
                repmat({'_'},length(seriesNums),1),...
                arrayfun(@num2str, numCells(seriesNums)', 'UniformOutput', 0)',...
                repmat({'_'},length(seriesNums),1),...
                cellType1Name(seriesNums)',...
                repmat({'_'},length(seriesNums),1),...
                cellType2Name(seriesNums)');
    end


    % Loop through each video within the dataset dataset

    for iSeries = seriesNums

        % Update a counter
            seriesCounter = seriesCounter+1;

        % Save the information
            CadCombo_BestFitExp(seriesCounter) = 1*contains(cellType1Name{iSeries},'Cdh3')+2*contains(cellType1Name{iSeries},'Cdh2');
            Day_BestFitExp(seriesCounter) = datasetNum;
            ReplicateWithinDay_BestFitExp(seriesCounter) = 1*(iSeries<=4)+2*(iSeries>4);
            ExpressionLevel_BestFitExp(seriesCounter) = 1*contains(cellType1Name{iSeries},'low')+2*contains(cellType1Name{iSeries},'high');

        % Select the path to the file
            expMatFilePath = [imageFilePath(1:(find(imageFilePath=='.',1,'last'))-1) ...
                sprintf('_S%1i_',iSeries) 'MaxIP_Segmentation.mat'];

        % Plot the comparison of the time lapse montage

            % Load the cell type segmentation data
                clear cellType_OverTime
                load(expMatFilePath,'cellType_OverTime');   

            % Pull out the file path for experimental image data
                prefixUS = find(expMatFilePath=='_',2,'last');
                prefixUS = prefixUS(1);
                expImagePath = [expMatFilePath(1:prefixUS) 'MaxIP_MATLAB.tiff'];     

            % Pull out the initial timepoint to analyze
                numTP = 1;

            % Read the file (but do not load the images)
                reader_MIP = bfGetReader(expImagePath);

            % Load the image for the experiments
                % Pull out the appropriate planes
                    iPlane1 = reader_MIP.getIndex(1 - 1, 1 -1, numTP - 1) + 1;
                    iPlane2 = reader_MIP.getIndex(1 - 1, 2 -1, numTP - 1) + 1;
                % Load the planes and convert to double precision values from 0 to 1
                    % where 1 maps to the largest possible for that specific
                    % data type (e.g., for unit16 1=2^16-1, for unit8 1=2^8-1)
                    cellType1ChannelIdxImage_MIP = im2double(bfGetPlane(reader_MIP, iPlane1));
                    cellType2ChannelIdxImage_MIP = im2double(bfGetPlane(reader_MIP, iPlane2));

            % Create a space to store the best fit values
                cellType_All{seriesCounter} = cellType_OverTime(:,:,1);
                rawImg_CellType1_All{seriesCounter} = cellType1ChannelIdxImage_MIP;
                rawImg_CellType2_All{seriesCounter} = cellType2ChannelIdxImage_MIP;

    end

end

%% Fig. 6 - Supp. Fig. 1 -- Plot the intensity

repNumOverall = nan(length(CadCombo_BestFitExp),1);
repNumOverall((Day_BestFitExp==1)&(ReplicateWithinDay_BestFitExp==1)) = 1;
repNumOverall((Day_BestFitExp==1)&(ReplicateWithinDay_BestFitExp==2)) = 2;
repNumOverall((Day_BestFitExp==2)&(ReplicateWithinDay_BestFitExp==1)) = 3;
repNumOverall((Day_BestFitExp==2)&(ReplicateWithinDay_BestFitExp==2)) = 4;

conditionNumOverall = nan(length(CadCombo_BestFitExp),1);
conditionNumOverall((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==1)) = 3;
conditionNumOverall((CadCombo_BestFitExp==1)&(ExpressionLevel_BestFitExp==2)) = 4;
conditionNumOverall((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==1)) = 1;
conditionNumOverall((CadCombo_BestFitExp==2)&(ExpressionLevel_BestFitExp==2)) = 2;

medianIntensity_CellType1 = nan(length(CadCombo_BestFitExp),1);
medianIntensity_CellType2 = nan(length(CadCombo_BestFitExp),1);


cmap = [1 0 0; 0 1 1];
markers = {'o','s','d','^'};


for seriesNum = 1:length(CadCombo_BestFitExp)

medianIntensity_CellType1(seriesNum) = median(rawImg_CellType1_All{seriesNum}(:));
medianIntensity_CellType2(seriesNum) = median(rawImg_CellType2_All{seriesNum}(:));

figure(8)
scatter(conditionNumOverall(seriesNum),medianIntensity_CellType1(seriesNum),...
    10,markers{ReplicateWithinDay_BestFitExp(seriesNum)+2*(Day_BestFitExp(seriesNum)-1)},'filled','MarkerFaceColor',...
    cmap(CadCombo_BestFitExp(seriesNum),:)./ExpressionLevel_BestFitExp(seriesNum), ...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
hold on;
scatter(conditionNumOverall(seriesNum)+0.3,medianIntensity_CellType2(seriesNum),...
    10,markers{ReplicateWithinDay_BestFitExp(seriesNum)+2*(Day_BestFitExp(seriesNum)-1)},'filled','MarkerFaceColor',...
    cmap(CadCombo_BestFitExp(seriesNum),:)./ExpressionLevel_BestFitExp(seriesNum), ...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerEdgeColor',[0 0 0])

end

figure(8)
hold off;
xlabel('Condition')
ylabel({'Median intensity','(first timepoint, unitless)'})
xticks([1 2 3 4])
xticklabels({'Cdh3-low+Cdh1_low','Cdh3-high+Cdh1-high','Cdh2-low+Cdh1-low','Cdh2-high+Cdh1-high'})
ylim([0 10^(-2)])
set(gca,'FontName','Helvetica','FontSize',5);
box on;

% Save a figure of the best fit
    % Pull out the file path
        figFilePath = [figOutputFolderPath 'Fig6SuppFig1_MeanIntensityByCondition'];
    % Choose the figure size and position
        figHandle = figure(8);
        figHandle.Position =  [488  438   150   135];
    % Save the file as a pdf
        exportgraphics(gcf,[figFilePath '.pdf'],'ContentType','vector')
    % Save as a png                
        saveas(gcf,[figFilePath '.png'])
