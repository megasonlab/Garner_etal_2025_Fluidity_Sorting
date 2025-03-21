
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
   
% Choose the path to the data (larger, lower-resolution scan)
    dataFolderPath = [simMasterFolderPath 'Fixed_Time\Results20240221\'];       

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
        load(matFilePath,'globalInfo','parameterValsNum','replicateNum')

    % Save the parameter values
        kT_All(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));
        kT_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,3);
        E_homo_All(matFileNum) = globalInfo.parameterVals1(globalInfo.paramCombos(parameterValsNum,1));
        E_homo_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,1);
        v_All(matFileNum) = globalInfo.parameterVals2(globalInfo.paramCombos(parameterValsNum,2));
        v_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,2);
        repNum_All(matFileNum) = replicateNum;


    end


end

%% Fig. 5a -- Plot the degree of sorting at 48hrs as a function of E_M and E_homo (showing tight coupling required for sorting)

    % Pull out the unique values for each parameter
        % Pull out the unique values for kT/E_homo
            unique_kTVals = unique(kT_All);
        % Pull out the unique values for kT/E_homo
            unique_E_homoVals = unique(E_homo_All);
        % Pull out the unique values for kT/E_homo
            unique_vVals = unique(v_All);

    % Choose the values for homotypic energy and viscosity to plot
        viscosityNum2Plot = 3;
        repNum2Plot = 1;
    % Pull out the subset of simulations to plot
        matfileNums2Plot = find((v_All==unique_vVals(viscosityNum2Plot))&(repNum_All==repNum2Plot));

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

        
    % Plot kT vs E_homo vs sorting
        figure(7)    
        scatter(kT_All(matfileNums2Plot).*E_homo_All(matfileNums2Plot),...
            E_homo_All(matfileNums2Plot),3,percentSorted_48hr_All(matfileNums2Plot'));
        hold on;
    % Plot the kT = E_homo line
       % plot(10.^([1.5:8.5]),10.^([1.5:8.5]),'k--')
        x2 = [10.^([1.5:8.5]), fliplr(10.^([1.5:8.5]))];
        inBetween = [10.^([1.5:8.5])*5, fliplr(10.^([1.5:8.5])/1.5)];
        fill(x2, inBetween, 'k','FaceAlpha',0.1,'EdgeAlpha',0);
        hold off;
    % Label the plot
        ylabel('E_{homo} (multiples of k_BT)')
        xlabel('E_M (multiples of k_BT)')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        colormap(cool*0.9)
        colorbar
        clim([0.5 1])
        xlim([10^1.5 10^8.5])
        ylim([10^1.5 10^8.5])
       % legend({'Sorting',sprintf('E_M = E_{homo}'),sprintf('Narrow window  \n where sorting can occur')},'Location','Southeast')
        legend({'Simulation data',sprintf('Narrow window \n around E_M = E_{homo} \n where sorting can occur')},'Location','Southeast')
        set(gca,'FontName','Helvetica','FontSize',5)


     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig5a_SortingEVsKTCoupling'];
        % Choose the figure size and position
            figHandle = figure(7);
            figHandle.Position =  [250   375   250   160];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

