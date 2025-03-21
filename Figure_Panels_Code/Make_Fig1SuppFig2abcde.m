
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

%% Plot the domain size fitting method

% Choose the path to the data
    simDataFolderPath = [simMasterFolderPath 'MC\Results20241017\'];

% Pull out the mat files
    matFiles = dir([simDataFolderPath '*_out.mat']);
    % Sort the files the way a human would
       [matFiles] = natsortfiles(matFiles);

% Load the number of timepoints
    load([simDataFolderPath matFiles(1).name],'numTP2Save');

% Choose the mat file number and timepoint to plot
    matFileNum = 10;
    numTP = round(numTP2Save/2);

% Determine the matFilePath
    matFilePath = [simDataFolderPath matFiles(matFileNum).name];

% Load the file
    clear cellType_onGrid_OverTime
    clear sizeX
    load(matFilePath,'cellType_onGrid_OverTime','sizeX')

% Calculate the domain size as the peak in radially-averaged 2D FFT of 2D spatial autocorrelation

    % Pull out the raw cell type data
        cellType = squeeze(cellType_onGrid_OverTime(:,:,numTP));
    % Calculate the 2D autocorrelation function
        cellType_2Dacorr = autocorr2d(cellType);
        separationInCells = -floor(sizeX/2):(floor(sizeX/2)-1);
    % Calculate the 2D FFT of the 2D autocorrelation function    
        cellType_2Dacorr_2DFFT = fft2(cellType_2Dacorr)/sizeX;
    % Calculate the radially-averaged 2D FFT of the 2D autocorrelation function 
        [cellType_2Dacorr_2DFFT_RMean,Y_STD_ThisTP] = calculateRadialAverage(...
            abs(cellType_2Dacorr_2DFFT),[1 1]);
        % Extract the single-sided version
            cellType_2Dacorr_2DFFT_RMean = cellType_2Dacorr_2DFFT_RMean(1:round((size(cellType,1)*2/3)+1));   
    % Extract the wavenumbers associated with the Fourier transform
        waveNumbers = 0:(size(cellType,1)*2/3);
    % Calculate the mean wavenumber
        power = 1;
        meanK_2DAcorr_2DFFT = sum(waveNumbers.*cellType_2Dacorr_2DFFT_RMean.^power')./sum(cellType_2Dacorr_2DFFT_RMean.^power);

    % Plot the steps
        figure(1)
        subplot(1,4,1)
        imagesc(cellType)
        xlabel('Position (# cells)')
        ylabel('Position (# cells)')
        colorbar
        title('Cell type (\phi(x,y))')
        subplot(1,4,2)
        imagesc(separationInCells,separationInCells,cellType_2Dacorr)
        xlabel('Separation (# cells)')
        ylabel('Separation (# cells)')
        colorbar
        title({'2D spatial autocorrelation function','C(r,t) = <\phi(x+r,t)\phi(x,t)>','(equal time pair correlation function)'})
        subplot(1,4,3)
        imagesc(abs(cellType_2Dacorr_2DFFT))
        xlabel('Separation (# cells)')
        ylabel('Separation (# cells)')
        colorbar
        title({'2D FFT of','2D spatial autocorrelation function','S(k,t) = FFT(C(r,t))','(equal time structure factor)'})
        xlabel('Wavenumber (k_x)')
        ylabel('Wavenumber (k_y)')
        colorbar
        subplot(1,4,4)
        scatter(waveNumbers, cellType_2Dacorr_2DFFT_RMean,30,'g','filled')
        hold on;
        yl = ylim();
        plot([meanK_2DAcorr_2DFFT meanK_2DAcorr_2DFFT],yl,'m--','LineWidth',2)
        hold off;
        legend({'',sprintf('Mean k = %2.2f',meanK_2DAcorr_2DFFT)})
        title({'Radially-averaged 2D FFT of','2D spatial autocorrelation function','S(k,t) = FFT(C(r,t))','(equal time structure factor)'})
        ylabel('Mean FFT')
        xlabel('Wavenumber (k)')
        drawnow

 % Save the figure 
    % Choose the filename and path for the figure
        destinationTrack = [figOutputFolderPath 'Fig1SuppFig2abcd_DomainSizeCalculation'];
    % Choose the figure size and position
        figHandle = figure(1);
        figHandle.Position =  [20 350 1500 300];
    % Save the file as a pdf
        exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
    % Save as a png                
        saveas(gcf,[destinationTrack '.png'])

%% Plot the three simulations for the lowest temperature on the same plot

% Choose the kT value to plot
    kTVal2Plot = 0.25;
 
% Load the global swap MC simululations
   
    % Choose the path to the data
        simDataFolderPath = [simMasterFolderPath 'MC\Results20241017\'];
    
    % Pull out the mat files
        matFiles = dir([simDataFolderPath '*_out.mat']);
        % Sort the files the way a human would
           [matFiles] = natsortfiles(matFiles);

    % Load the number of timepoints
        clear numTP2Save
        load([simDataFolderPath matFiles(1).name],'numTP2Save')

    % For recording the parameter inputs
        kT_All_MCGlobal = nan([1 length(matFiles)]); 
    % Time points
        timeInSteps_OverTime_MCGlobal = nan([length(matFiles),numTP2Save]);
    % Pre-allocate space to store the mean domain size over time   
        meanK_2DAcorr_2DFFT_All_MCGlobal = nan([length(matFiles),numTP2Save]);

    % Find the apprpriate simulation and fit the domain size over time 
    for matFileNum = 1:length(matFiles)
    
        % Pull out the file path
            matFilePath = [simDataFolderPath matFiles(matFileNum).name];
    
        % Load the file
            clear globalInfo
            clear parameterValsNum
            load(matFilePath,'globalInfo','parameterValsNum')
        
        % Save the parameter values
            kT_All_MCGlobal(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));

        if any(kT_All_MCGlobal(matFileNum)==kTVal2Plot)

            % Load the rest of the tile
                clear timeInSteps_OverTime
                clear cellType_onGrid_OverTime
                clear sizeX
                load(matFilePath,'timeInSteps_OverTime','cellType_onGrid_OverTime','sizeX')

            % Choose the timpoints to analyze
                timepoints2Analyze = 1:30:find(~isnan(cellType_onGrid_OverTime(1,1,:)),1,'last');
    
            % Choose whether to plot the fit
                doPlot=false;
    
            % Run the analysis
                [meanK_2DAcorr_2DFFT_All_MCGlobal(matFileNum,timepoints2Analyze)] = ...
                    calculateDomainSizeFromStructureFactor(cellType_onGrid_OverTime,timepoints2Analyze,sizeX,doPlot);
            % Save the time information
                timeInSteps_OverTime_MCGlobal(matFileNum,timepoints2Analyze) = timeInSteps_OverTime(timepoints2Analyze)+1;
            % Save the sysetm size
                sizeX_MCGlobal = sizeX;
       
        end
    end

% Load the local swap MC simululations
   
    % Choose the path to the data
        simDataFolderPath = [simMasterFolderPath 'MC\Results20241017_1\'];
    
    % Pull out the mat files
        matFiles = dir([simDataFolderPath '*_out.mat']);
        % Sort the files the way a human would
           [matFiles] = natsortfiles(matFiles); 

    % Load the number of timepoints
        clear numTP2Save
        load([simDataFolderPath matFiles(1).name],'numTP2Save')

    % For recording the parameter inputs
        kT_All_MCLocal = nan([1 length(matFiles)]); 
    % Time points
        timeInSteps_OverTime_MCLocal = nan([length(matFiles),numTP2Save]);
    % Pre-allocate space to store the mean domain size over time   
        meanK_2DAcorr_2DFFT_All_MCLocal = nan([length(matFiles),numTP2Save]);

    % Find the apprpriate simulation and fit the domain size over time 
    for matFileNum = 1:length(matFiles)
    
        % Pull out the file path
            matFilePath = [simDataFolderPath matFiles(matFileNum).name];
    
        % Load the file
            clear globalInfo
            clear parameterValsNum
            load(matFilePath,'globalInfo','parameterValsNum')

        % Save the parameter values
            kT_All_MCLocal(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));

        if any(kT_All_MCLocal(matFileNum)==kTVal2Plot)

            % Load the rest of the tile
                clear timeInSteps_OverTime
                clear cellType_onGrid_OverTime
                clear sizeX
                load(matFilePath,'timeInSteps_OverTime','cellType_onGrid_OverTime','sizeX')

            % Choose the timpoints to analyze
                timepoints2Analyze = 1:30:find(~isnan(cellType_onGrid_OverTime(1,1,:)),1,'last');
    
            % Choose whether to plot the fit
                doPlot=false;
    
            % Run the analysis
                [meanK_2DAcorr_2DFFT_All_MCLocal(matFileNum,timepoints2Analyze)] = ...
                    calculateDomainSizeFromStructureFactor(cellType_onGrid_OverTime,timepoints2Analyze,sizeX,doPlot);
            % Save the time information
                timeInSteps_OverTime_MCLocal(matFileNum,timepoints2Analyze) = timeInSteps_OverTime(timepoints2Analyze)+1;
            % Save the sysetm size
                sizeX_MCLocal = sizeX;


          
        end
    end

% Load the fixed time step simululations
   
    % Choose the path to the data
        simDataFolderPath = [simMasterFolderPath 'Fixed_Time\Results20240819_1\'];  
    
    % Pull out the mat files
        matFiles = dir([simDataFolderPath '*_out.mat']);
        % Sort the files the way a human would
           [matFiles] = natsortfiles(matFiles); 

    % Load the number of timepoints
        clear numTP2Save
        load([simDataFolderPath matFiles(1).name],'numTP2Save')

    % For recording the parameter inputs
        kT_All_DT = nan([1 length(matFiles)]); 
    % Time points
        timeInSteps_OverTime_DT = nan([length(matFiles),numTP2Save]);
    % Pre-allocate space to store the mean domain size over time   
        meanK_2DAcorr_2DFFT_All_DT = nan([length(matFiles),numTP2Save]);

    % Find the apprpriate simulation and fit the domain size over time 
    for matFileNum = 1:length(matFiles)
    
        % Pull out the file path
            matFilePath = [simDataFolderPath matFiles(matFileNum).name];
    
        % Load the file
            clear globalInfo
            clear parameterValsNum
            load(matFilePath,'globalInfo','parameterValsNum')

        % Save the parameter values
            kT_All_DT(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));

        if any(kT_All_DT(matFileNum)==kTVal2Plot)

            % Load the rest of the tile
                clear timeInMin_OverTime
                clear cellType_onGrid_OverTime
                clear sizeX
                load(matFilePath,'timeInMin_OverTime','cellType_onGrid_OverTime','sizeX')

            % Choose the timpoints to analyze
                timepoints2Analyze = 1:30:find(~isnan(cellType_onGrid_OverTime(1,1,:)),1,'last');
    
            % Choose whether to plot the fit
                doPlot=false;
    
            % Run the analysis
                [meanK_2DAcorr_2DFFT_All_DT(matFileNum,timepoints2Analyze)] = ...
                    calculateDomainSizeFromStructureFactor(cellType_onGrid_OverTime,timepoints2Analyze,sizeX,doPlot);
            % Save the time information
                timeInSteps_OverTime_DT(matFileNum,timepoints2Analyze) = timeInMin_OverTime(timepoints2Analyze);
            % Save the sysetm size
                sizeX_DT = sizeX;

      
        end
    end


% Plot the result    
    figure(2)
    % Global swaps
        scatter(mean(timeInSteps_OverTime_MCGlobal(kT_All_MCGlobal==kTVal2Plot,timepoints2Analyze),1,'omitnan'),...
            sizeX./mean(meanK_2DAcorr_2DFFT_All_MCGlobal(kT_All_MCGlobal==kTVal2Plot,timepoints2Analyze),1,'omitnan'),30,'b*')
        hold on;
    % Local swaps
        scatter(mean(timeInSteps_OverTime_MCLocal(kT_All_MCLocal==kTVal2Plot,timepoints2Analyze),1,'omitnan'),...
            sizeX./mean(meanK_2DAcorr_2DFFT_All_MCLocal(kT_All_MCLocal==kTVal2Plot,timepoints2Analyze),1,'omitnan'),30,'bo')
    % Discrete time
        scatter(mean(timeInSteps_OverTime_DT(kT_All_DT==kTVal2Plot,:),1,'omitnan')*60000,...
            sizeX./mean(meanK_2DAcorr_2DFFT_All_DT(kT_All_DT==kTVal2Plot,:),1,'omitnan'),30,'bd')
    % Analytical theory
        timeTheory = 10.^(2.5:0.1:5);
        plot(timeTheory,0.7*10^(0)*timeTheory.^(1/3),'k--')
    % Clean up the plot
        hold off;
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        xlabel('Time (steps for MC, ms for DT)')
        ylabel('Domain size (# cells, 1/<k>)')
        set(gca,'ytick',[1 3 10 30 100])
        set(gca,'xtick',10.^(1:1:6))
        ylim([2.5 sizeX_DT])
        xlim([10^(1) 10^6])
        legend({'MC (global)';'MC (local)';'DT (local)';'N^{1 / 3}'},'Location','Northwest')

 % Save the figure 
    % Choose the filename and path for the figure
        destinationTrack = [figOutputFolderPath 'Fig1SuppFig2e_GlobalvsLocalvsDT'];
    % Choose the figure size and position
        figHandle = figure(2);
        figHandle.Position =  [500 500 350 275];
    % Save the file as a pdf
        exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
    % Save as a png                
        saveas(gcf,[destinationTrack '.png'])