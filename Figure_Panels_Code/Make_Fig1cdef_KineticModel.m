
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

%% Fig. 1c -- Plot the probability of being sorted as a function of the  
% adhesion energy (E_A) for different choices of the motility energy (E_M)

% Choose the set of parameters to plot
    % Choose the adhesion energies for the given state (in multiples of
    % k_BT_L)
        E_kT_L_Vals = [-2:0.01:0];
    % Choose the motility energies (in multiples of
    % k_BT_L)
        kT_M_kT_L_Vals = [0.1 0.3 1 3 100]; 
    % Choose the thermal energy of the lab temperature in pN nm
        kT_L_pNnm = 4.2;
    % Choose the viscosity in pN ms / nm
        v_pNms_nm_Vals = 1000*10^(-3);    
    % Choose the distance between cell centers in nm
        L_nm = 10*10^3;
    % Choose the colormap
        CM = jet(length(kT_M_kT_L_Vals))*0.8;
        CM(4,:) = CM(5,:);
        CM(5,:) = [1 0 0]*0.8;
  
    % Pull up the figure
        figure(1)
    % Plot the sorting probabilty vs differential adhesion energy for
    % different choices of kT
        for kT_Num = 1:length(kT_M_kT_L_Vals)    
        plot(abs(E_kT_L_Vals), (60*10^3)*(kT_M_kT_L_Vals(kT_Num)*kT_L_pNnm./(v_pNms_nm_Vals*L_nm^2)).*...
            exp(E_kT_L_Vals./kT_M_kT_L_Vals(kT_Num)),'--','LineWidth',2,...
                'Color',[CM(kT_Num,:) 1])
        hold on;
        end
        hold off;
    % Label the plot
        set(gca,'YScale','log')
      %  ylim([0 1])
        xlabel('Total adhesion energy (E_A)')
        ylabel('Fluidity (neighbor exchange rate, min^{-1})')
        title('Fluidity vs E_M and E_A')
        legend([cellstr(num2str(kT_M_kT_L_Vals', 'E_M=%1.1f'))],...
            'Location','Southwest')
        set(gca,'FontName','Helvetica','FontSize',5); 

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig1c_FluidityVskTAndEhomo'];
        % Choose the figure size and position
            figHandle = figure(1);
            figHandle.Position =  [475   325   210   160];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

%% Load the simulation dataset

% Extract the data folder path
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

%% Fig. 1d -- Plot a time lapse montage of sorting for an example simulation

% Choose the example simulation to plot

    % Calculate the relevant ratios
        ERatio_All = kT_All;
        VRatio_All = E_homo_All./v_All;

    % Choose the target parameter values
        ERatio_Target = 0.85;
        VRatio_Target = 2*10^5;
    
    % Find the closest parameter values in the dataset
        [ERatio_Target_Closest_Diff, ERatio_Target_Closest_Idx] = min(abs(ERatio_All-ERatio_Target)',[],1);
        ERatio_Target_Closest = ERatio_All(ERatio_Target_Closest_Idx);
        [VRatio_Target_Closest_Diff, VRatio_Target_Closest_Idx] = min(abs(VRatio_All-VRatio_Target)',[],1);
        VRatio_Target_Closest = VRatio_All(VRatio_Target_Closest_Idx);
    
    % Pull out the mat file number for the simuations of interest
        matFileNums2Plot = find((ERatio_All==ERatio_Target_Closest)&(VRatio_All==VRatio_Target_Closest));
        matFilIdx2Plot = 4;
        matFileNum2Plot = matFileNums2Plot(matFilIdx2Plot);
    
    % Pull out the file path
        matFilePath = [simDataFolderPath matFiles(matFileNum2Plot).name];

% Plot the time lapse montage

    % Load the relevant data
        clear cellID_onGrid_OverTime
        clear cellType_onGrid_OverTime
        clear timeInMin_OverTime
        load(matFilePath,'cellID_onGrid_OverTime','cellType_onGrid_OverTime','timeInMin_OverTime')
    
    % Choose the time points to plot in min
        timeInMin2Plot_Ideal = [0 10 3*60 (60*48)];
    
    % Find the associated time in the simulations that are closest to the
    % timepoints selected for plotting
        % Calculate the time different between all possible pairs of
        % the experimental time and the actual timepoints
            timeDifferenceAllPairs = abs(timeInMin_OverTime(:)-timeInMin2Plot_Ideal(:).');
        % Find the closest simulated timepoint to each experimental
        % timepoint
            [M, timePointIndex] = min(timeDifferenceAllPairs,[],1);
        % Pull out the actual time values
            timeInMin2Plot = timeInMin_OverTime(timePointIndex);
    
    % Plot the time lapse montage
    for timePointNum = 1:length(timeInMin2Plot)
        figure(2)
        subplot(2,length(timeInMin2Plot),timePointNum)
        imagesc(cellID_onGrid_OverTime(:,:,timePointIndex(timePointNum)))
        % Title the panel with the time
            if timeInMin2Plot(timePointNum)==0
            title(sprintf('Time = %1.1f min',timeInMin2Plot(timePointNum))) 
            elseif timeInMin2Plot(timePointNum)<0.4
            title(sprintf('Time = %1.1f sec',timeInMin2Plot(timePointNum)*60)) 
            elseif timeInMin2Plot(timePointNum)<59
            title(sprintf('Time = %1.1f min',timeInMin2Plot(timePointNum))) 
            else
            title(sprintf('Time = %1.1f hrs',timeInMin2Plot(timePointNum)/60)) 
            end
        xlabel('Number of cells - X')
        ylabel('Number of cells - Y')
        set(gca,'FontName','Helvetica','FontSize',5);
        subplot(2,length(timeInMin2Plot),length(timeInMin2Plot) + timePointNum)
        imagesc(cellType_onGrid_OverTime(:,:,timePointIndex(timePointNum)))
        xlabel('Number of cells - X')
        ylabel('Number of cells - Y')
        set(gca,'FontName','Helvetica','FontSize',5);
    end
    
     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig1d_SortingTimeLapseMontage'];
        % Choose the figure size and position
            figHandle = figure(2);
            figHandle.Position =  [250   375   415   180];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

%% Fig. 1e -- Plot the associated sorting dynamics

    for matFileNum2Plot = matFileNums2Plot(matFilIdx2Plot)
    
    % Pull out the file path
        matFilePath = [simDataFolderPath matFiles(matFileNum2Plot).name];

    % Load the relevant data
        clear numSameCellTime_OverTime
        clear numNeighbors_ByPosition
        clear timeInMin_OverTime
        clear domainWidthinMicrons_RepByPixel
        load(matFilePath,'numSameCellTime_OverTime','numNeighbors_ByPosition','timeInMin_OverTime','domainWidthinMicrons_RepByPixel')

    % Plot the degree of sorting over time

        % Calculate the percent sorted
            fraction_Sorted = nanmean(numSameCellTime_OverTime./numNeighbors_ByPosition,1);
    
        % Plot the degree of sorting over time
            figure(3)
            yyaxis left
            scatter(timeInMin_OverTime,fraction_Sorted*100,1,'filled')
            hold on;
            yyaxis right
            scatter(timeInMin_OverTime,domainWidthinMicrons_RepByPixel.Mean./domainWidthinMicrons_RepByPixel.cellDiameterInUM,1,'filled')
            hold on;
        
    end
    % Label the plot
        yyaxis left
        hold off
        yyaxis right
        hold off
       % set(gca,'XScale','log')
       % set(gca,'YScale','log')
        xlabel('Time (min)')
        yyaxis left
        ylabel({'Percent sorted','(% same cell type neighbors)'})
        ylim([50 100])
        yyaxis right
        ylabel({'Degree of Sorting','(Domain size, # cells)'})
        ylim([1.5 5])
        xlim([10^(-2) 60*49])
        set(gca,'FontName','Helvetica','FontSize',5)

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig1e_SortingVsTime'];
        % Choose the figure size and position
            figHandle = figure(3);
            figHandle.Position =  [250   375   150   100];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

%% Fig. 1f -- Plot the associated fluidity dynamics

    % Load the relevant data
        clear sizeGrid
        clear cellRadius
        clear numCells
        clear numTP2Save
        clear cellID_onGrid_OverTime
        clear timeInMin_OverTime
        load(matFilePath,'sizeGrid','cellRadius','numCells','numTP2Save','cellID_onGrid_OverTime','timeInMin_OverTime')

    % Convert the data to standard format
        
        % Reshape the data to be numTracks by XYZ
        % Preallocate space to store the data
            trackPositionsXYZ_UM_reshaped = nan([numCells,3,numTP2Save]);
        
        % Determine the X,Y,Z positions of each position on the grid
            [xPosInUM,yPosInUM] = meshgrid((1:sizeGrid(1))*cellRadius/1000,(1:sizeGrid(2))*cellRadius/1000);
        
        % Lopo through each timepoint and record the cell positions
            for tpNum = 1:numTP2Save
                cellID_Temp = cellID_onGrid_OverTime(:,:,tpNum);
                trackPositionsXYZ_UM_reshaped(cellID_Temp(:),1,tpNum) = ...
                    xPosInUM(:);
                trackPositionsXYZ_UM_reshaped(cellID_Temp(:),2,tpNum) = ...
                    yPosInUM(:);
            end
            % Record 0 for z
                trackPositionsXYZ_UM_reshaped(:,3,:) = 0;
        
        % Create the time vector
            time_Min_reshaped = repmat(timeInMin_OverTime', [numCells,1]);

    % Run the neighbor exchange analysis
        [neighborsExchangedThisTP, slopes, time2Fit, window, slopeMedian, time_Min_NeighborExchange] = ...
            calculateNeighborExchangeRate(trackPositionsXYZ_UM_reshaped,time_Min_reshaped,cellRadius/1000*1.1*ones([numTP2Save 1]));
    
    % Calculate the summary statistics
        % Calculate the cumulative number of neighbors exchanged
            cumNumNeighborsExchanged = cumsum(neighborsExchangedThisTP,2);
        % Calculate the cumulative number of neighbors exchanged
            meanCumNumNeighborsExchanged = mean(cumsum(neighborsExchangedThisTP,2),1,'omitnan');
        % Calculate the time step
            timeStep_min = diff(time_Min_reshaped,1,2);
        % Calculate the neighbor exchange rate
            neighborExchangeRate = neighborsExchangedThisTP./timeStep_min;
        % Calculate the mean number of neighbors exchanged between each timepoint
            meanNumNeighborsExchangedThisTP = median(neighborsExchangedThisTP,1,'omitnan');
        % Calculate the mean neigbor exchange rate in each timepoint
            meanNeighborExchangeRateThisTP = meanNumNeighborsExchangedThisTP./mean(timeStep_min,1);

    % Remove the small values
        slopesNonZero = slopes;
        slopesNonZero(slopesNonZero<10^(-10)) = nan;

    % Plot the cumulative number of neighbors exchanged over time
        figure(4)
        subplot(1,2,1)
        % Choose the data color
            alphaVal = 1/sqrt(sqrt(numCells));
            experimentColor = [0.1 0.1 0.9 alphaVal]; 
        % Plot all of the curves simulteaneously
            plot(time_Min_NeighborExchange,cumNumNeighborsExchanged,'color',experimentColor)
            hold on;
        % Plot the mean
            plot(time_Min_NeighborExchange,meanCumNumNeighborsExchanged,'-o',...
                'LineWidth',1,'MarkerSize',1,'color',[0 1 0]*0.7)
            hold off;
        % Label the plot
            ylabel('Cumulative neighbors exchanged')
            xlabel('Time (min)')
            title('Cumulative neighbor exchange over time')
            xlim([min(time_Min_NeighborExchange) max(time_Min_NeighborExchange)])
            set(gca,'FontName','Helvetica','FontSize',5)

    % Plot the neighbor exchange rate
        figure(4)
        subplot(1,2,2)
    % Plot the slope of the median
        plot(time_Min_NeighborExchange(time2Fit),slopeMedian,'-o',...
            'LineWidth',1,'MarkerSize',1,'color',[0 1 0]*0.7)
    % Label the plot
        hold off;
        ylabel({'Neighbor exchange rate','(neighbors per min)'})
        xlabel('Time (min)')
        title({'Moving linear regression of','cumulative neighbor exchange'}) 
        legend(cellstr(num2str(window, 'Window=%-d timepoints')))
        set(gca,'FontName','Helvetica','FontSize',5)

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig1f_FluidityVsTime'];
        % Choose the figure size and position
            figHandle = figure(4);
            figHandle.Position =  [250   375   350   135];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

