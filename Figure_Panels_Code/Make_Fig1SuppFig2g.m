
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

%% Plot the steady state simulations with 50 and 150 box size on the same plot

% Choose the kT value to plot
    kTVal2Plot = 1.3;
 
% Load the larger box size simulations
   
    % Choose the path to the data
        simDataFolderPath = [simMasterFolderPath 'Steady_State\Results20240927\'];

    % Pull out the mat files
        matFiles = dir([simDataFolderPath '*_out.mat']);
        % Sort the files the way a human would
           [matFiles] = natsortfiles(matFiles); 

    % Load the number of timepoints
        clear numTP2Save
        load([simDataFolderPath matFiles(1).name],'numTP2Save')

    % For recording the parameter inputs
        kT_All_Large = nan([1 length(matFiles)]); 
    % Time points
        timeInSteps_OverTime_Large = nan([length(matFiles),numTP2Save]);
    % Pre-allocate space to store the mean domain size over time   
        meanK_2DAcorr_2DFFT_All_Large = nan([length(matFiles),numTP2Save]);

    % Find the apprpriate simulation and fit the domain size over time 
    for matFileNum = 1:length(matFiles)
    
        % Pull out the file path
            matFilePath = [simDataFolderPath matFiles(matFileNum).name];
    
        % Load the file
            clear globalInfo
            clear parameterValsNum
            load(matFilePath,'globalInfo','parameterValsNum')

        % Save the parameter values
            kT_All_Large(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));

        if any(abs(kT_All_Large(matFileNum)-kTVal2Plot)<0.01)

            % Load the rest of the tile
                clear timeInMin_OverTime
                clear cellType_onGrid_OverTime
                clear sizeX
                load(matFilePath,'timeInMin_OverTime','cellType_onGrid_OverTime','sizeX','BCs')
                BCs

            % Choose the timpoints to analyze
                timepoints2Analyze = 1:30:find(~isnan(cellType_onGrid_OverTime(1,1,:)),1,'last');
    
            % Choose whether to plot the fit
                doPlot=false;
    
            % Run the analysis
                [meanK_2DAcorr_2DFFT_All_Large(matFileNum,timepoints2Analyze)] = ...
                    calculateDomainSizeFromStructureFactor(cellType_onGrid_OverTime,timepoints2Analyze,sizeX,doPlot);
            % Save the time information
                timeInSteps_OverTime_Large(matFileNum,timepoints2Analyze) = timeInMin_OverTime(timepoints2Analyze);
            % Save the sysetm size
                sizeX_Large = sizeX;

      
        end
    end

% Load the smaller box size simulations
   
    % Choose the path to the data
        simDataFolderPath = [simMasterFolderPath 'Steady_State\Results20240221\'];
    
    % Pull out the mat files
        matFiles = dir([simDataFolderPath '*_out.mat']);
        % Sort the files the way a human would
           [matFiles] = natsortfiles(matFiles); 

    % Load the number of timepoints
        clear numTP2Save
        load([simDataFolderPath matFiles(1).name],'numTP2Save')

    % For recording the parameter inputs
        kT_All_Small = nan([1 length(matFiles)]); 
    % Time points
        timeInSteps_OverTime_Small = nan([length(matFiles),numTP2Save]);
    % Pre-allocate space to store the mean domain size over time   
        meanK_2DAcorr_2DFFT_All_Small = nan([length(matFiles),numTP2Save]);

    % Find the apprpriate simulation and fit the domain size over time 
    for matFileNum = 1:length(matFiles)
    
        % Pull out the file path
            matFilePath = [simDataFolderPath matFiles(matFileNum).name];
    
        % Load the file
            clear globalInfo
            clear parameterValsNum
            load(matFilePath,'globalInfo','parameterValsNum')

        % Save the parameter values
            kT_All_Small(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));

        if any(abs(kT_All_Small(matFileNum)-kTVal2Plot)<0.01)

            % Load the rest of the tile
                clear timeInMin_OverTime
                clear cellType_onGrid_OverTime
                clear sizeX
                load(matFilePath,'timeInMin_OverTime','cellType_onGrid_OverTime','sizeX','BCs')
                BCs

            % Choose the timpoints to analyze
                timepoints2Analyze = 1:30:find(~isnan(cellType_onGrid_OverTime(1,1,:)),1,'last');
    
            % Choose whether to plot the fit
                doPlot=false;
    
            % Run the analysis
                [meanK_2DAcorr_2DFFT_All_Small(matFileNum,timepoints2Analyze)] = ...
                    calculateDomainSizeFromStructureFactor(cellType_onGrid_OverTime,timepoints2Analyze,sizeX,doPlot);
            % Save the time information
                timeInSteps_OverTime_Small(matFileNum,timepoints2Analyze) = timeInMin_OverTime(timepoints2Analyze);
            % Save the sysetm size
                sizeX_Small = sizeX;

      
        end
    end

%% Plot the result    
    figure(1)
    % Small box size
        scatter(mean(timeInSteps_OverTime_Small(kT_All_Small==kTVal2Plot,:),1,'omitnan'),...
            sizeX_Small./mean(meanK_2DAcorr_2DFFT_All_Small(kT_All_Small==kTVal2Plot,:),1,'omitnan'),30,'filled','mo')
        hold on;
    % Large box size
        scatter(mean(timeInSteps_OverTime_Large(kT_All_Large==kTVal2Plot,:),1,'omitnan'),...
            sizeX_Large./mean(meanK_2DAcorr_2DFFT_All_Large(kT_All_Large==kTVal2Plot,:),1,'omitnan'),30,'filled','go')
    % Clean up the plot
        hold off;
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        xlabel('Time (min)')
        ylabel('Domain size (# cells, 1/<k>)')
        set(gca,'ytick',[1 3 10 30 100])
        set(gca,'xtick',10.^(-6:2:6))
        ylim([2.5 sizeX_Large])
        xlim([10^(-4) 10^4])
        legend({sprintf('Grid Size = %i',sizeX_Small);sprintf('Grid Size = %i',sizeX_Large)},'Location','Northwest')

 % Save the figure 
    % Choose the filename and path for the figure
        destinationTrack = [figOutputFolderPath 'Fig1SuppFig2g_GridSize'];
    % Choose the figure size and position
        figHandle = figure(1);
        figHandle.Position =  [500 450 275 250];
    % Save the file as a pdf
        exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
    % Save as a png                
        saveas(gcf,[destinationTrack '.png'])