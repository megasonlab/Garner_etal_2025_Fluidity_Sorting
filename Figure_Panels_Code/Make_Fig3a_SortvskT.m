
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

%% Load an example dataset

% Choose the path to the data (vary kT_eff)
    dataFolderPath = [simMasterFolderPath 'Fixed_Time\Results20240819_1\'];      

% Pull out the mat files
    matFiles = dir([dataFolderPath '*_out.mat']);
    % Sort the files the way a human would
       [matFiles] = natsortfiles(matFiles);

% Load the number of timepoints
    load([dataFolderPath matFiles(1).name],'numTP2Save');

% For recording the parameter inputs
    kT_All = nan([1 length(matFiles)]);  
    kT_All_Idx = nan([1 length(matFiles)]);  
    E_homo_All = nan([1 length(matFiles)]);  
    E_homo_All_Idx = nan([1 length(matFiles)]);  
    v_All = nan([1 length(matFiles)]);  
    v_All_Idx = nan([1 length(matFiles)]);  
    repNum_All = nan([1 length(matFiles)]); 
    rSwap_0_All = nan([1 length(matFiles)]); 

% Loop through each file and fit the domain size over time and the domain
% growth law
for matFileNum = 1:length(matFiles)

    % Print the file number
      %  matFileNum

    % Pull out the file path
        matFilePath = [dataFolderPath matFiles(matFileNum).name];

    try

    % Load the file
        clear globalInfo
        clear parameterValsNum
        load(matFilePath,'globalInfo','parameterValsNum','replicateNum','rSwap_0')

    % Save the parameter values
        kT_All(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));
        kT_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,3);
        E_homo_All(matFileNum) = globalInfo.parameterVals1(globalInfo.paramCombos(parameterValsNum,1));
        E_homo_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,1);
        v_All(matFileNum) = globalInfo.parameterVals2(globalInfo.paramCombos(parameterValsNum,2));
        v_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,2);
        repNum_All(matFileNum) = replicateNum;
        rSwap_0_All(matFileNum) = rSwap_0;


    end


end

% Calculate the rate of swapping for each pair (min^(-1))
    meanSwapRatePerMin_All = rSwap_0_All.*60.*...
        exp(-4./kT_All);

%% Plot sorting vs time for representative values of kT  and endpoint sorting vs kT

% Choose which parameters to plot

    % Pull out the unique values for each parameter
        % Pull out the unique values for kT/E_homo
            unique_kTVals = unique(kT_All);
        % Pull out the unique values for kT/E_homo
            unique_E_homoVals = unique(E_homo_All);
        % Pull out the unique values for kT/E_homo
            unique_vVals = unique(v_All);
    
    % Choose the values for homotypic energy and viscosity to plot
        EhomoNum2Plot = 1;
        viscosityNum2Plot = 1;
        repNum2Plot = 1;
    % Choose the representative values of kT to emphasize
        kTVals2Emphasize = [0.25 0.5 0.75 1 1.25 1.5]; 
    
    % Pull out the subset of simulations to plot
        matfileNums2Plot = find((E_homo_All==unique_E_homoVals(EhomoNum2Plot))&...
            (v_All==unique_vVals(viscosityNum2Plot))&(repNum_All==repNum2Plot));

% For these files, calculate the % sorted
    % Pre-allocate space to store the domain size at each timepoint
        timeInMin_All = nan([length(matFiles) numTP2Save]);
        percentSorted_All = nan([length(matFiles) numTP2Save]);
    % Loop through each file and fit the domain size over time and the domain
    % growth law
    for matFileNum = matfileNums2Plot
        % Pull out the file path
            matFilePath = [dataFolderPath matFiles(matFileNum).name];
        try
        % Load the relevant data
            clear numSameCellTime_OverTime
            clear numNeighbors_ByPosition
            load(matFilePath,'numSameCellTime_OverTime','numNeighbors_ByPosition','timeInMin_OverTime')    
        % Calculate the percent sorted
            percentSorted_All(matFileNum,:) = nanmean(numSameCellTime_OverTime./numNeighbors_ByPosition,1);
            timeInMin_All(matFileNum,:) = timeInMin_OverTime;
        end  
    end    
    % Calculate the degree of sorting at the 48hr timepoint for each parameter
    % set
        percentSorted_48hr_All = nanmean(percentSorted_All(:,(end-50):end),2);

% Plot the the probabiltiy of sorting vs kT for the chosen values for
% E_homo and viscosity
    figure(5)
    subplot(1,3,1)
    scatter(kT_All(matfileNums2Plot),...%.*E_homo_All(matfileNums2Plot),...
        meanSwapRatePerMin_All(matfileNums2Plot),3,...
        ...%kT_All(matfileNums2Plot),...%.*E_homo_All(matfileNums2Plot),...
        log10(meanSwapRatePerMin_All(matfileNums2Plot)),...
        'filled')
    hold on;
    subplot(1,3,2)
    scatter(kT_All(matfileNums2Plot),...%.*E_homo_All(matfileNums2Plot),...
        percentSorted_48hr_All(matfileNums2Plot)*100,3,...
        ...%kT_All(matfileNums2Plot),...%.*E_homo_All(matfileNums2Plot),...
        log10(meanSwapRatePerMin_All(matfileNums2Plot)),...
        'filled')
    hold on;


    for kTNum2Emphasize = 1:length(kTVals2Emphasize)
        
        % Choose the subset to plot (the closest value for kT)
            matfileNum2Plot = find((E_homo_All==unique_E_homoVals(EhomoNum2Plot))&...
                (v_All==unique_vVals(viscosityNum2Plot))&...
                (abs(kT_All-kTVals2Emphasize(kTNum2Emphasize))==min(abs(kT_All-kTVals2Emphasize(kTNum2Emphasize))))&...
                (repNum_All==repNum2Plot));


        % Add star to fluidity vs kT graph
            subplot(1,3,1)
            scatter(kT_All(matfileNum2Plot),...%.*E_homo_All(matfileNum2Plot),...
                meanSwapRatePerMin_All(matfileNum2Plot),10,...
                ...%kT_All(matfileNum2Plot),...%.*E_homo_All(matfileNum2Plot),...
               log10(meanSwapRatePerMin_All(matfileNum2Plot)),...
               'pentagram','filled','MarkerEdgeColor',[0 0 0])
        % Add star to sorting vs kT graph
            subplot(1,3,2)
            scatter(kT_All(matfileNum2Plot),...%.*E_homo_All(matfileNum2Plot),...
                percentSorted_48hr_All(matfileNum2Plot)*100,10,...
                ...%kT_All(matfileNum2Plot),...%.*E_homo_All(matfileNum2Plot),...
                log10(meanSwapRatePerMin_All(matfileNum2Plot)),...
                'pentagram','filled','MarkerEdgeColor',[0 0 0])
        % Add plot of sorting vs time
            subplot(1,3,3)
            scatter(timeInMin_All(matfileNum2Plot,:),...
                percentSorted_All(matfileNum2Plot,:)*100,3,...
                ...%kT_All(matfileNum2Plot).*...%E_homo_All(matfileNum2Plot)*...
                log10(meanSwapRatePerMin_All(matfileNum2Plot)).*...
                ones(size(timeInMin_All(matfileNum2Plot,:))),...
                'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
                hold on;

    end

% Choose the color limits
    climVals = [-3.5 2.25];

% Label the plots
    % Fluidity 
        subplot(1,3,1)
        hold off;
        box on;
        % Color
            colormap(jet*0.9)
            clim(climVals)
        % Y-axis
            ylim([10^(-7.5) 10^2])
            set(gca,'YTick',[10^(-7) 10^(-5) 10^(-3) 10^(-1) 10^1 10^3 10^5])
            set(gca,'YScale','log')
            ylabel({'Neighbor exchange rate' 'at 50% sorted (min^{-1})'})
        % X-axis
            set(gca,'XTick',[0 0.5 1 1.5 2])
            xlim([-0.1 2.1])
            xlabel('E_M (multiples of E_{homo})')
        % Parameter labels
            xL=xlim;
            yL=ylim;
            text(0.99*xL(2),0.99*yL(2),sprintf('E_{homo} = %1.1g k_BT_L \n viscosity = %1.1d Pa s',...
                unique_E_homoVals(EhomoNum2Plot),unique_vVals(viscosityNum2Plot)),'HorizontalAlignment','right','VerticalAlignment','top',...
                'FontName','Helvetica','FontSize',5)
        % Font
            set(gca,'FontName','Helvetica','FontSize',5)
    % Sorting 
        subplot(1,3,2)
        hold off;
        box on;
        % Color
            colormap(jet*0.9)
            clim(climVals)
        % Y-axis
            set(gca,'YTick',[50 60 70 80 90 100])
            ylim([49 100])
            ylabel({'Sorting at 48 hrs','(% same cell type)'})
        % X-axis
            set(gca,'XTick',[0 0.5 1 1.5 2])
            xlim([-0.1 2.1])
            xlabel('E_M (multiples of E_{homo})')
        % Paarameter labels
            xL=xlim;
            yL=ylim;
            text(0.99*xL(2),0.99*yL(2),sprintf('E_{homo} = %1.1g k_BT_L \n viscosity = %1.1d Pa s',...
                unique_E_homoVals(EhomoNum2Plot),unique_vVals(viscosityNum2Plot)),'HorizontalAlignment','right','VerticalAlignment','top',...
                'FontName','Helvetica','FontSize',5)
            set(gca,'FontName','Helvetica','FontSize',5)
    % Sorting over time
        subplot(1,3,3)
        hold off;
        box on;
        % Color
            colormap(jet*0.9)
            clim(climVals)
        % Y-axis
            ylim([49 100])
            set(gca,'YTick',[50 60 70 80 90 100])
            ylabel({'Sorting','(% same cell type)'})
        % X-axis
            xlim([10^(-4) 25*60])
            set(gca,'XScale','log')    
            xlabel('Time (min)')
        % Parameter labels
            xL=xlim;
            yL=ylim;
            text(0.99*xL(2),0.99*yL(2),sprintf('E_{homo} = %1.1g k_BT_L \n viscosity = %1.1d Pa s',...
                unique_E_homoVals(EhomoNum2Plot),unique_vVals(viscosityNum2Plot)),'HorizontalAlignment','right','VerticalAlignment','top',...
                'FontName','Helvetica','FontSize',5)
        % Text
            set(gca,'FontName','Helvetica','FontSize',5)

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig3a_SortingVisualVskT'];
        % Choose the figure size and position
            figHandle = figure(5);
            figHandle.Position =  [250   375   375   85];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

