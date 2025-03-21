
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

%% Load an the steady state dataset
   
% Choose the path to the data (steady state)
    dataFolderPath = [simMasterFolderPath 'Steady_State\Results20240221\'];  

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

% For recording raw data
    % Time points
        timeInMin_All = cell(length(matFiles),1);
    % Fraction of neighbors sorted
        fracSameCellType_All = cell(length(matFiles),1);

% For recording the fit results
    % Fit object
        fitResults_All = cell([1 length(matFiles)]); 
    % Fit parameters
        fracSameCellTypeAtEq = nan([1 length(matFiles)]);
        charTimescaleInMin = nan([1 length(matFiles)]); 
        hillCoeff = nan([1 length(matFiles)]); 
        asymmetryCoeff = nan([1 length(matFiles)]); 
    % Time to steady state
        time2SteadyStateInMin = nan([1 length(matFiles)]);

% Record how long it takes to load the data
tic

% Loop through each file and fit the domain size over time and the domain
% growth law
for matFileNum = 1:length(matFiles)

    % Print the file number
      %  matFileNum

    % Pull out the file path
        matFilePath = [dataFolderPath matFiles(matFileNum).name];

   % try

    % Load the file
        clear globalInfo
        clear parameterValsNum
        clear numSameCellTime_OverTime
        clear numNeighbors_ByPosition
        clear timeInMin_OverTime
        load(matFilePath,'globalInfo','parameterValsNum','replicateNum',...
            'numSameCellTime_OverTime','numNeighbors_ByPosition','timeInMin_OverTime')

    % Save the parameter values
        kT_All(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));
        kT_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,3);
        E_homo_All(matFileNum) = globalInfo.parameterVals1(globalInfo.paramCombos(parameterValsNum,1));
        E_homo_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,1);
        v_All(matFileNum) = globalInfo.parameterVals2(globalInfo.paramCombos(parameterValsNum,2));
        v_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,2);
        repNum_All(matFileNum) = replicateNum;


     % Save the raw data
        timeInMin_All{matFileNum} = timeInMin_OverTime;
        fracSameCellType_All{matFileNum} = mean(numSameCellTime_OverTime./...
                numNeighbors_ByPosition,1)';

    % Perform the fit
        t2Fit = timeInMin_All{matFileNum};
        sort2Fit =  fracSameCellType_All{matFileNum};
        idx2fit = ~isnan(sort2Fit);
        [fitresult] = createFitAsymHill(t2Fit(idx2fit),sort2Fit(idx2fit),false);

    % Record the fit results
        % Save the whole thin
            fitResults_All{matFileNum} = fitresult; 
        % Pull out the best fit parameters
            bestFitParams = coeffvalues(fitresult);
        % Separate them into the relevant variables
            fracSameCellTypeAtEq(matFileNum) = bestFitParams(1);
            charTimescaleInMin(matFileNum) = bestFitParams(2);
            hillCoeff(matFileNum) = bestFitParams(3);
            asymmetryCoeff(matFileNum) = bestFitParams(4);

    % Calculate the time to reach steady state    
        % The system asymptotically approches equillibrium. 
        % Choose how close the system has to be to the equillibrium value 
        % to call it equillibrium
            percentSortedWindow = 0.001;
        % Calculate the time to steady state
            [time2SteadyStateInMin(matFileNum)] = ...
                evalInverseAsymHill(fracSameCellTypeAtEq(matFileNum),...
                charTimescaleInMin(matFileNum),hillCoeff(matFileNum),...
                asymmetryCoeff(matFileNum), ...
                fracSameCellTypeAtEq(matFileNum)-percentSortedWindow);

end

%%

% Choose which parameters to plot

    % Pull out the unique values for each parameter
        % Pull out the unique values for kT/E_homo
            unique_kTVals = unique(kT_All);
        % Pull out the unique values for kT/E_homo
            unique_E_homoVals = unique(E_homo_All);
        % Pull out the unique values for kT/E_homo
            unique_vVals = unique(v_All);
    % Choose twhic parameters to plot
        kTNum2Plot = 3;
        repNum2Plot = 1;

% Plot a single curve and its fit
    figure(8)
    subplot(1,4,1)

    % Pull out the subset of simulations to plot
        matfileNums2Plot = find((kT_All==unique_kTVals(kTNum2Plot))&(repNum_All==repNum2Plot));

    % Plot the results
        % Choose the timepoints to plot the prediction
            time2plotPredict = 10.^(linspace(-6,log10(2*time2SteadyStateInMin(matfileNums2Plot)),1000));
        % Calculate the best fit line in this region
            percentSortedPredict = feval(fitResults_All{matfileNums2Plot},time2plotPredict);
        % Plot the raw data
            scatter(timeInMin_All{matfileNums2Plot},fracSameCellType_All{matfileNums2Plot},3,...
                unique_kTVals(kTNum2Plot)*ones(size(timeInMin_All{matfileNums2Plot})),...
                'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
            hold on;
        % Plot the best fit line
            p1 = plot(time2plotPredict,percentSortedPredict,'k--',...
                'LineWidth',1);
            set(p1,'HandleVisibility','off') % Don't put it on the legend
        % Plot the steady state sorting time
            set(gca,'XScale','log')
            xl = xlim();
            yl = ylim();
            plot(time2SteadyStateInMin(matfileNums2Plot)*[1 1],yl,'r-.')
            plot(xl, fracSameCellTypeAtEq(matfileNums2Plot)*[1 1],'r--')
            hold off;
        % Clean up the plot
            colormap(jet*0.9)
            clim([0 1.6])%*10^log10(unique_E_homoVals(EhomoNum2Plot)))
            xlim([time2plotPredict(1) time2plotPredict(end)])
            ylim([0.5 1]) 
            xlabel('Time (min)')
            ylabel({'Degree of sorting at 48 hrs','(% same cell type neighbors)'})
            legend({'Simulation data','Best fit to sigmoid'},'Location','Northwest')
            set(gca,'FontName','Helvetica','FontSize',5)

% Plot the time to steady state

    % Choose the sorting times to highlight
        % times2Highlight = [1 60 (60*24) (60*24*7) (60*24*30) (60*24*365) (60*24*365*10) (60*24*365*100) (60*24*365*1000)];
        % labels = {'1 min','1 hr','1 day','1 week','1 month','1 year','10 years','100 years','1000 years'};
        times2Highlight = [1 60 (60*24) (60*24*30) (60*24*365*10) (60*24*365*1000)];
        labels = {'1 min','1 hr','1 day','1 month','10 years','1000 years'};
    
    % Plot the sorting rate vs the motility energy
        % Plot the sorting
            figure(8)
            subplot(1,4,2)
            scatter(kT_All,1./time2SteadyStateInMin,3,kT_All,'filled')
            hold on;
        % Add labels for specific timepoints
            predictedKTAtHighlight=[];
            [kT_All_Sorted I_sorted] = sort(kT_All); 
            time2SteadyStateInMin_Sorted = log(time2SteadyStateInMin(I_sorted));
            for textLabelNum = 1:length(times2Highlight)
                [closeTimeVals closeIdxs] = mink(abs(log(times2Highlight(textLabelNum))-time2SteadyStateInMin_Sorted)',3);
                predictedKTAtHighlight(textLabelNum) = mean(kT_All_Sorted(closeIdxs));
            end
            plot(predictedKTAtHighlight,1./times2Highlight,'k*')
            text(predictedKTAtHighlight+0.06,1./times2Highlight,labels,...
                'VerticalAlignment','top','HorizontalAlignment','left',...
                'FontWeight','bold','FontSize',5,'BackgroundColor', [0.9 0.9 0.9 0.51])
            hold off;
    
        % Clean up the plot
            set(gca,'YScale','log')
            colormap(jet*0.9)
            clim([0 1.6])%*10^log10(unique_E_homoVals(EhomoNum2Plot)))
            xlim([0.9 2.1])
            ylim([10^(-11) 10^(2)])
            ylabel({'Sorting rate (min^{-1})'})
            xlabel({'E_M/E_{homo}'})
            set(gca,'FontName','Helvetica','FontSize',5)

% Plot the degree of sorting at steady state
    
    % Plot the sorting rate vs the motility energy
        % Plot the sorting
            figure(8)
            subplot(1,4,3)
            scatter(kT_All,fracSameCellTypeAtEq,3,kT_All,'filled')
            hold on;
    
        % Clean up the plot
            colormap(jet*0.9)
            clim([0 1.6])%*10^log10(unique_E_homoVals(EhomoNum2Plot)))
            xlim([0.9 2.1])
            ylim([0.5 1])
            ylabel({'Degree of sorting at 48 hrs','(% same cell type neighbors)'})
            xlabel({'E_M/E_{homo}'})
            set(gca,'FontName','Helvetica','FontSize',5)

% Plot the degree of sorting at steady state vs the sorting time
    
    % Plot the sorting rate vs the motility energy
        % Plot the sorting
            figure(8)
            subplot(1,4,4)
            scatter(1./time2SteadyStateInMin,fracSameCellTypeAtEq,3,kT_All,'filled')
            hold on;
    
        % Clean up the plot
            set(gca,'XScale','log')
            colormap(jet*0.9)
            clim([0 1.6])%*10^log10(unique_E_homoVals(EhomoNum2Plot)))
            xlim([10^(-11) 10^(2)])
            ylim([0.5 1])
            ylabel({'Degree of sorting at 48 hrs','(% same cell type neighbors)'})
            xlabel({'Sorting rate (min^{-1})'})
            set(gca,'FontName','Helvetica','FontSize',5)

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig3bcde_SteadyStateSortingVsEM'];
        % Choose the figure size and position
            figHandle = figure(8);
            figHandle.Position =  [250   375   765   135];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

%% Plot many examples of the best fit 

    % Prepare the figure
        figure(88)
        numCols = round(sqrt(length(unique_kTVals)));
        numRows = ceil(length(unique_kTVals)/numCols);

    % Loop through all the kT vals and show the fit
    for kTNum2Plot = 1:length(unique_kTVals)
        subplot(numRows,numCols,kTNum2Plot)

    % Pull out the subset of simulations to plot
        matfileNums2Plot = find((kT_All==unique_kTVals(kTNum2Plot))&(repNum_All==repNum2Plot));

    % Plot the results
        % Choose the timepoints to plot the prediction
            time2plotPredict = 10.^(linspace(-6,log10(2*time2SteadyStateInMin(matfileNums2Plot)),1000));
        % Calculate the best fit line in this region
            percentSortedPredict = feval(fitResults_All{matfileNums2Plot},time2plotPredict);
        % Plot the raw data
            scatter(timeInMin_All{matfileNums2Plot},fracSameCellType_All{matfileNums2Plot}*100,3,...
                unique_kTVals(kTNum2Plot)*ones(size(timeInMin_All{matfileNums2Plot})),...
                'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
            hold on;
        % Plot the best fit line
            p1 = plot(time2plotPredict,percentSortedPredict*100,'k--',...
                'LineWidth',1);
        % Plot the steady state sorting time
            set(gca,'XScale','log')
            xl = xlim();
            yl = ylim();
            hold off;
        % Clean up the plot
            colormap(jet*0.9)
            clim([0 1.6])%*10^log10(unique_E_homoVals(EhomoNum2Plot)))
            xlim([time2plotPredict(1) time2plotPredict(end)])
            ylim([49 100]) 
            xlabel('Time (min)')
            ylabel({'Sorting','(% same cell type)'})
            set(gca,'ytick',[50 60 70 80 90 100])
            minXTick = log10(time2plotPredict(1));
            maxXTick = ceil(log10(time2plotPredict(end)));
            set(gca,'xtick',10.^(minXTick:ceil((maxXTick-minXTick)/4):maxXTick))
            title(sprintf('E_M = %0.2g E_{homo}',unique_kTVals(kTNum2Plot)))
            set(gca,'FontName','Helvetica','FontSize',5)
    end

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig3SuppFig1_PlotFitsForAllCurves'];
        % Choose the figure size and position
            figHandle = figure(88);
            figHandle.Position =  [50   115   700   555];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

   