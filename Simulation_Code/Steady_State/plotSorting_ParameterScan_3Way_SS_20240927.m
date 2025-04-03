
% Clear the system
    close all;
    clear all;

% Choose the dataset to plot    
    % Choose the path to the data
        dataFolderPath = ['C:\Users\Rikki\Documents\Postdoc\Research\'...
            'Modeling\Sorting_Fluidity\Parameter_Scans\Steady_State\Results20240927\'];       
    % Choose the template for sprintf
        templateFor2Sprintf = 'kT = %1.2f';

% Load the relevant information
    % Pull out the mat files
        matFiles = dir([dataFolderPath '*_out.mat']);
      %  matFiles = dir([dataFolderPath '*32248720_*_out.mat']);
    % Sort the files the way a human would
       [matFiles] = natsortfiles(matFiles);

% Create variables to store the data
    % For recording raw data
        % Time points
            timeInMin_All = cell(length(matFiles),1);
        % Fraction of neighbors sorted
            fracSameCellType_All = cell(length(matFiles),1);
        % Create a legend text variable
            legendText = cell(length(matFiles),1);
    % For recording the fit results
        % Fit object
            fitResults_All = cell([1 length(matFiles)]); 
        % Fit parameters
            fracSameCellTypeAtEq = nan([1 length(matFiles)]);
            charTimescaleInMin = nan([1 length(matFiles)]); 
            hillCoeff = nan([1 length(matFiles)]); 
            asymmetryCoeff = nan([1 length(matFiles)]); 
    % For recording the fit results
        sortState_48hr = nan([1 length(matFiles)]); 
        kT_All = nan([1 length(matFiles)]);  
        kT_All_Idx = nan([1 length(matFiles)]);  

% Initialize a counter
    legendTextCounter=1;


% Pull up the figure
    f = figure(200);

% Loop through each file and plot the data alongside the best fit
for matFileNum = 1:length(matFiles)

    %try

    % Load the file
        clear globalInfo
        clear parameterValsNum
        clear numTP2Save
        clear timeInMin_OverTime
        clear numSameCellTime_OverTime
        clear numNeighbors_ByPosition
        load([dataFolderPath matFiles(matFileNum).name],...
            'globalInfo','parameterValsNum','numTP2Save',...
            'timeInMin_OverTime','numSameCellTime_OverTime',...
            'numNeighbors_ByPosition')


    % Close any graphics objects
   %     all_figs = findobj(0, 'type', 'figure');
    %    delete(setdiff(all_figs, 200));

    % Choose the timepoints to plot
        TPs2Plot = 1:numTP2Save;
        timeInMin = timeInMin_OverTime(TPs2Plot);

    % Adjust whichever parameters we're scanning through
        % Homotypic energy (in # kT_lab)
            E_homo = -globalInfo.parameterVals1(globalInfo.paramCombos(parameterValsNum,1));
        % Effective kT as a fraction of E_homo
            kT_eff = -globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3))*E_homo;
        % Choose the viscosity of the tissue in pN nm ms
            etaVal = globalInfo.parameterVals2(globalInfo.paramCombos(parameterValsNum,2))*10^(-3);
    
    % Count the number of conditions
        % Columns represent a different values for E_homo
            numCols = length(globalInfo.parameterVals1);
        % Rows represent a different values for viscosity
            numRows = length(globalInfo.parameterVals2);
        % Color represents different values of kT
        %    cmap = jet(length(globalInfo.parameterVals3))*0.9;
            cmapKTVals = linspace(0,2,100);
            cmap = jet(length(cmapKTVals))*0.9;

            % Calculate the time different between all possible pairs
                kTDifferenceAllPairs = abs(cmapKTVals'...
                    -globalInfo.parameterVals3);
            % Find the cloests value in the colormap
                [M, I] = min(kTDifferenceAllPairs,[],1);
            % Truncate the colormap
                cmap = cmap(I,:);
    % 
    % % Update the legend information and counter
    %    legendText{legendTextCounter} = sprintf(templateFor2Sprintf,...
    %        globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3)));
    %    legendTextCounter = legendTextCounter+1;

     % Save the raw data
        timeInMin_All{matFileNum} = timeInMin;
        fracSameCellType_All{matFileNum} = mean(numSameCellTime_OverTime(:,TPs2Plot)./...
                numNeighbors_ByPosition,1)';

    % Save the % sorted at the end of teh sim
        sortState_48hr(matFileNum) = nanmean(fracSameCellType_All{matFileNum}((end-100):end));
        kT_All(matFileNum) = globalInfo.parameterVals3(globalInfo.paramCombos(parameterValsNum,3));
        kT_All_Idx(matFileNum) = globalInfo.paramCombos(parameterValsNum,3);

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

    % Plot the results
        % Choose the timepoints to plot the prediction
            time2plotPredict = 10.^(linspace(-6,15,1000));
        % Calculate the best fit line in this region
            percentSortedPredict = feval(fitresult,time2plotPredict);
        % Plot the raw data
            scatter(timeInMin_All{matFileNum},fracSameCellType_All{matFileNum},...
                'filled','MarkerFaceColor',cmap(globalInfo.paramCombos(parameterValsNum,3),:),'MarkerFaceAlpha',0.2)
            hold on;
        % Plot the best fit line
            p1 = plot(time2plotPredict,percentSortedPredict,'-',...
                'LineWidth',1,'Color',cmap(globalInfo.paramCombos(parameterValsNum,3),:));
            set(p1,'HandleVisibility','off') % Don't put it on the legend
        % Clean up the plot
            set(gca,'XScale','log')
        %    set(gca,'YScale','log')
            ylim([0.5 1]) 
            % title(sprintf('E_{homo}=%1.1d',...
            %     globalInfo.parameterVals1(globalInfo.paramCombos(parameterValsNum,1))))
            % ylabel(sprintf('Viscosity=%1.1d',globalInfo.parameterVals2(globalInfo.paramCombos(parameterValsNum,2))))
   %         xlim([0 time2plotPredict(end)])

   % end

end

% Clean up the plot
figure(200);
%legend(legendText(~cellfun('isempty',legendText)))
hold off;
drawnow
xlim([10^(-7) 10^(7)])
ylabel('Mean fraction of same cell type neighbors')
xlabel('Time (min)')
title('Sorting dynamics as a function of kT_{eff}')      
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');

    % % Choose the video file path
    %     figFilePath = [dataFolderPath 'SortingDynamics.png'];
    % 
    %     % Save the figure
    %         % Choose the figure size and position
    %             figHandle = figure(200);
    %             figHandle.Position =  [500   250   550   425];
    %             set(figHandle,'PaperPositionMode','auto')
    %             print(figHandle,figFilePath,'-dpng',['-r' num2str(150)])

%% Calculate the time to reach equillibrium

    % The system asymptotically approches equillibrium. 
    % Choose how close the system has to be to the equillibrium value 
    % to call it equillibrium
        percentSortedWindow = 0.001;

    % Pre-allocate space to store the value
        time2SteadyStateInMin = nan([1 length(matFiles)]);

    % Loop through each file and calculate the time to equillibrium
    for matFileNum = 1:length(matFiles)
    
        [time2SteadyStateInMin(matFileNum)] = ...
            evalInverseAsymHill(fracSameCellTypeAtEq(matFileNum),...
            charTimescaleInMin(matFileNum),hillCoeff(matFileNum),...
            asymmetryCoeff(matFileNum), ...
            fracSameCellTypeAtEq(matFileNum)-percentSortedWindow);
    end


% Sort the kT values
    [kT_All_Sorted,I_sorted] = sort(kT_All);

    kT_All_Idx_Sorted = kT_All_Idx(I_sorted);

    good2plot = (1:length(I_sorted));
    [minVal minIdx] = min(time2SteadyStateInMin);
    good2plot(I_sorted==minIdx) = [];

% Choose the sorting times to highlight
times2Highlight = [1 60 (60*24) (60*24*7) (60*24*30)];% (60*24*365) (60*24*365*10) (60*24*365*100)];
labels = {'1 min','1 hr','1 day','1 week','1 month',};%'1 year','10 years','100 years'};

figure(10)
subplot(2,2,1)
scatter(kT_All_Sorted(good2plot),fracSameCellTypeAtEq(I_sorted(good2plot)),...
    20, cmap(kT_All_Idx_Sorted(good2plot),:),'filled')
ylabel({'Fraction of same cell type neighbors','at equillibrium'})
xlabel({'kT_{eff}/E_{homo}'})
ylim([0.5 1])
xlim([1 2])
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');

subplot(2,2,2)
scatter(kT_All_Sorted(good2plot),1./time2SteadyStateInMin(I_sorted(good2plot)),...
    30, cmap(kT_All_Idx_Sorted(good2plot),:), 'filled')
hold on;
predictedKTAtHighlight=[];
timeVals = log(time2SteadyStateInMin(I_sorted(good2plot)));
kTValsTemp = kT_All_Sorted(good2plot);
for textLabelNum = 1:length(times2Highlight)
    [closeTimeVals closeIdxs] = mink(abs(log(times2Highlight(textLabelNum))-timeVals)',3);
    predictedKTAtHighlight(textLabelNum) = mean(kTValsTemp(closeIdxs));
end
plot(predictedKTAtHighlight,1./times2Highlight,'k*')
text(predictedKTAtHighlight-0.06,1./times2Highlight,labels,...
    'VerticalAlignment','bottom','HorizontalAlignment','right')
hold off;
ylabel({'Sorting rate (min^{-1})'})
xlabel({'kT_{eff}/E_{homo}'})
%ylim([0.5 1])
xlim([1 2.1])
ylim([10^(-6) 2*10^(1)])
set(gca,'YScale','log')
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');

colororder([0.5 0 0.5;0 0.5 0.5])
subplot(2,2,3)
yyaxis left
scatter(kT_All_Sorted(good2plot),1./time2SteadyStateInMin(I_sorted(good2plot)),...
    30, cmap(kT_All_Idx_Sorted(good2plot),:), 'filled',...
    'MarkerEdgeColor',[0.4 0.1 0.5],'LineWidth',1)
    %'LineWidth',1.5) 
hold on;
predictedKTAtHighlight=[];
timeVals = log(time2SteadyStateInMin(I_sorted(good2plot)));
kTValsTemp = kT_All_Sorted(good2plot);
for textLabelNum = 1:length(times2Highlight)
    [closeTimeVals closeIdxs] = mink(abs(log(times2Highlight(textLabelNum))-timeVals)',3);
    predictedKTAtHighlight(textLabelNum) = mean(kTValsTemp(closeIdxs));
end
plot(predictedKTAtHighlight,1./times2Highlight,'k*')
text(predictedKTAtHighlight-0.06,1./times2Highlight,labels,...
    'VerticalAlignment','bottom','HorizontalAlignment','right',...
    'FontWeight','bold','FontSize',12,'BackgroundColor', [0.9 0.9 0.9 0.51])
hold off;
set(gca, 'YColor', [0.4 0.1 0.5]);
ylabel({'Sorting rate (\bullet min^{-1})'})
xlabel({'kT_{eff}/E_{homo}'})
%ylim([0.5 1])
xlim([1 2.1])
ylim([10^(-6) 2*10^(1)])
set(gca,'YScale','log')
yyaxis right
scatter(kT_All_Sorted(good2plot),fracSameCellTypeAtEq(I_sorted(good2plot)),...
    30, cmap(kT_All_Idx_Sorted(good2plot),:),'^','filled',...
    'MarkerEdgeColor',[0 0.4 0.8],'LineWidth',0.25)
    %'LineWidth',1.25)
set(gca, 'YColor', [0 0.4 0.8]);
ylabel({'Fraction of same cell type neighbors','at equillibrium (\Delta)'})
xlabel({'kT_{eff}/E_{homo}'})
ylim([0.5 1])
xlim([0.95 2.05])
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');


subplot(2,2,4)
scatter(fracSameCellTypeAtEq(I_sorted(good2plot)),1./time2SteadyStateInMin(I_sorted(good2plot)),...
    20, cmap(kT_All_Idx_Sorted(good2plot),:), 'filled')
hold on;
predictedSortingAtHighlight=[];
timeVals = log(time2SteadyStateInMin(I_sorted(good2plot)));
sortVals = fracSameCellTypeAtEq(I_sorted(good2plot));
for textLabelNum = 1:length(times2Highlight)
    [closeTimeVals closeIdxs] = mink(abs(log(times2Highlight(textLabelNum))-timeVals)',3);
    predictedSortingAtHighlight(textLabelNum) = mean(sortVals(closeIdxs));
end
plot(predictedSortingAtHighlight,1./times2Highlight,'k*')
text(predictedSortingAtHighlight,1./times2Highlight,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
hold off;
xlabel({'Fraction of same cell type neighbors','at equillibrium'})
ylabel({'Sorting rate (min^{-1})'})
set(gca,'YScale','log')
xlim([0.5 1])
ylim([10^(-6) 2*10^(1)])
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');

% Choose the video file path
    figFilePath = [dataFolderPath 'SortingTradeoff.png'];

    % Save the figure
        % Choose the figure size and position
            figHandle = figure(10);
            figHandle.Position =  [350   1000   1400   1000];
            set(figHandle,'PaperPositionMode','auto')
            print(figHandle,figFilePath,'-dpng',['-r' num2str(150)])

%%

figure(300)
subplot(2,2,1)
scatter(kT_All_Sorted(good2plot),fracSameCellTypeAtEq(I_sorted(good2plot)),...
    20, cmap(kT_All_Idx_Sorted(good2plot),:), 'filled')
ylabel({'Fraction of same cell type neighbors','at equillibrium'})
xlabel({'kT_{eff}/E_{homo}'})
ylim([0.5 1])
xlim([1 2])
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
subplot(2,2,2)
scatter(kT_All_Sorted(good2plot),charTimescaleInMin(I_sorted(good2plot)),...
    20, cmap(kT_All_Idx_Sorted(good2plot),:), 'filled')
ylabel({'Characteristic timescale (min)'})
xlabel({'kT_{eff}/E_{homo}'})
%ylim([0.5 1])
xlim([1 2])
set(gca,'YScale','log')
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
subplot(2,2,3)
scatter(kT_All_Sorted(good2plot),hillCoeff(I_sorted(good2plot)),...
    20, cmap(kT_All_Idx_Sorted(good2plot),:), 'filled')
ylabel({'Steepness (unitless)','(Hill coefficient)'})
xlabel({'kT_{eff}/E_{homo}'})
%ylim([0.5 1])
xlim([1 2])
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
subplot(2,2,4)
scatter(kT_All_Sorted(good2plot),asymmetryCoeff(I_sorted(good2plot)),...
    20, cmap(kT_All_Idx_Sorted(good2plot),:), 'filled')
ylabel({'Asymmetry'})
xlabel({'kT_{eff}/E_{homo}'})
%ylim([0.5 1])
xlim([1 2])
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');

% 
%     sortState_48hr_Sorted = sortState_48hr(I_sorted);
%     kT_All_Idx_Sorted = kT_All_Idx(I_sorted);
% 
%     globalInfo.paramCombos(parameterValsNum,3)
% 
% 
% figure(201)
% %plot(kT_All_Sorted,sortState_48hr_Sorted,'.-')
% scatter(kT_All_Sorted,sortState_48hr_Sorted, 20, cmap(kT_All_Idx_Sorted,:), 'filled')
% xlabel('kT (relative to E_{homo})')
% ylabel('Mean fraction of same cell type neighbors')
% title('Sorting at 48 hrs')
% ylim([0.5 1])
% xlim([0 2.1])
% set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');

% 
% %% Plot the best fit data
% 
% close(figure(1))
% figure(1)
% subplot(2,3,1)
% plot(paramVals(1:length(matFiles)),fracSameCellTypeAtEq,...
%     'k-o','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% xlabel('kT (relative to E_{homo})')
% ylabel({'Fraction of same cell type neighbors','at equillibrium'})
% subplot(2,3,2)
% plot(paramVals(1:length(matFiles)),charTimescaleInMin,...
%     'k-o','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% xlabel('kT (relative to E_{homo})')
% ylabel('Characteristic timescale (min)')
% subplot(2,3,4)
% plot(paramVals(1:length(matFiles)),hillCoeff,...
%     'k-o','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% xlabel('kT (relative to E_{homo})')
% ylabel('Steepness')
% subplot(2,3,5)
% plot(paramVals(1:length(matFiles)),asymmetryCoeff,...
%     'k-o','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% xlabel('kT (relative to E_{homo})')
% ylabel('Asymmetry')
% 
% % Calculate the time to reach equillibrium
% 
%     % The system asymptotically approches equillibrium. 
%     % Choose how close the system has to be to the equillibrium value 
%     % to call it equillibrium
%         percentSortedWindow = 0.01;
% 
%     % Pre-allocate space to store the value
%         time2SteadyStateInMin = nan([1 length(matFiles)]);
% 
%     % Loop through each file and calculate the time to equillibrium
%     for matFileNum = 1:length(matFiles)
% 
%         [time2SteadyStateInMin(matFileNum)] = ...
%             evalInverseAsymHill(fracSameCellTypeAtEq(matFileNum),...
%             charTimescaleInMin(matFileNum),hillCoeff(matFileNum),...
%             asymmetryCoeff(matFileNum), ...
%             fracSameCellTypeAtEq(matFileNum)-percentSortedWindow);
%     end
% 
% figure(1)
% subplot(2,3,3)
% plot(paramVals(1:length(matFiles)),time2SteadyStateInMin,...
%     'k-o','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% xlabel('kT (relative to E_{homo})')
% ylabel({'Time to equillibirum (min)'})
% set(gca,'YScale','log')
% 
% 
% figure(1)
% subplot(2,3,6)
% colororder([0.65 0.5 0.8; 0 0.75 0.75])
% yyaxis left
% p1 = plot(paramVals(1:length(matFiles)),1./time2SteadyStateInMin,...
%     '-o','LineWidth',1','MarkerSize',2);
% % hold on;
% % plot(paramVals(1:length(matFiles)),1/60*ones([1 length(matFiles)]),'k-')
% % plot(paramVals(1:length(matFiles)),1/(60*24)*ones([1 length(matFiles)]),'k--')
% % plot(paramVals(1:length(matFiles)),1/(60*24*7)*ones([1 length(matFiles)]),'k.-')
% % plot(paramVals(1:length(matFiles)),1/(60*24*365*100)*ones([1 length(matFiles)]),'k*-')
% % hold off;
% set(p1, 'MarkerFaceColor', get(p1,'Color')); 
% xlabel('kT (relative to E_{homo})')
% ylabel({'Sorting rate (min^{-1})'})
% set(gca,'YScale','log')
% yyaxis right
% p2 = plot(paramVals(1:length(matFiles)),fracSameCellTypeAtEq,...
%     '-o','LineWidth',1','MarkerSize',2);
% set(p2, 'MarkerFaceColor', get(p2,'Color')); 
% xlabel('kT (relative to E_{homo})')
% ylabel({'Fraction of same cell type neighbors','at equillibrium'})
% 
% %%
% figure(2)
% plot(time2SteadyStateInMin,fracSameCellTypeAtEq,...
%     'k-o','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% hold on;
% times2Highlight = [60 (60*24) (60*24*7) (60*24*365) (60*24*365*10) (60*24*365*100)];
% labels = {'1 hr','1 day','1 week','1 year','10 years','100 years'};
% predictedSortingAtHighlight = interp1(log(time2SteadyStateInMin),fracSameCellTypeAtEq,log(times2Highlight));
% plot(times2Highlight,predictedSortingAtHighlight,'*')
% text(times2Highlight,predictedSortingAtHighlight,labels,'VerticalAlignment','top','HorizontalAlignment','left')
% hold off;
% ylabel({'Fraction of same cell type neighbors','at equillibrium'})
% xlabel({'Time to equillibirum (min)'})
% set(gca,'XScale','log')
% 
% %% Fit the first n timepoints for the last one
% 
% % Choose the file to plot
%     numMatFile2Plot = 5;
% 
% % Pull out the data for this file
%     timeVals = timeInMin_All{numMatFile2Plot};
%     sortVals = fracSameCellType_All{numMatFile2Plot};    
%     time2SteadyStateInMinVal = time2SteadyStateInMin(numMatFile2Plot);
% 
% % Choose have many timepoints to plot
%     numFits = 100;
% % Choose which endpoints to sample
%     endpoints2sample = 11:max(1,round(length(timeVals)/numFits)):length(timeVals);
%     numFits = length(endpoints2sample);
% % Create a colormap
%     CM = jet(numFits);
% 
% % Create a space to store the data
%     SS_Sorted_Predict = nan([numFits,1]);
%     Time_Sorted_Predict = nan([numFits,1]);
%     SS_Bool = false([numFits,1]);
% 
% 
% for numFit = 1:numFits
% 
%     % Pull out the timepoints to fit
%         tp2Fit = 1:endpoints2sample(numFit);
% 
%     % Perform the fit
%         [fitresult] = createFitAsymHill(timeVals(tp2Fit), sortVals(tp2Fit), false);
% 
%     % Pull out the best fit parameters
%         bestFitParams = coeffvalues(fitresult);
%     % Separate them into the relevant variables
%         SS_Sorted_Predict(numFit) = bestFitParams(1);
%         Time_Sorted_Predict(numFit) = ...
%             evalInverseAsymHill(bestFitParams(1),...
%             bestFitParams(2),bestFitParams(3),...
%             bestFitParams(4), ...
%             bestFitParams(1)-percentSortedWindow);
% 
%     % % Plot the results
%     %     % Pull up the figure
%     %     figure(400);
%     %     % numRows = floor(sqrt(numFits));
%     %     % numCols = ceil(numFits/numRows);
%     %   %  subplot(numRows,numCols,numFit)
%     %     % Choose the timepoints to plot the prediction
%     %         time2plotPredict = 10.^(linspace(-6,15,1000));
%     %     % Calculate the best fit line in this region
%     %         percentSortedPredict = feval(fitresult,time2plotPredict);
%     %     % Plot the raw data
%     %     if numFit==1
%     %         scatter(timeVals,sortVals,...
%     %             'filled','MarkerFaceColor','k','MarkerFaceAlpha',0.2)
%     %     end
%     %         hold on;
%     %     % Plot the best fit line
%     %         p1 = plot(time2plotPredict,percentSortedPredict,'-','LineWidth',1,'Color',CM(numFit,:));
%     %         set(p1,'HandleVisibility','off') % Don't put it on the legend
%     %     % Clean up the plot
%     %         set(gca,'XScale','log')
%     %         set(gca,'YScale','log')
%     %         ylim([0.5 1]) 
%     %         xlim([0 time2plotPredict(end)])
% 
%     % Figure out whether the simulation has reached steady state
%         SS_Bool(numFit) = (nanmean(sortVals(tp2Fit((end-10):end)))>= (bestFitParams(1)-0.01));
% 
% 
% 
% end
% figure(400)
% hold off;
% ylabel('Mean fraction of same cell type neighbors')
% xlabel('Time (min)')
% 
% %%
% 
% 
%       SS_Time_Bool = (timeVals(endpoints2sample)'>(10*Time_Sorted_Predict));
% 
% figure(500)
% plot(timeVals(endpoints2sample),SS_Sorted_Predict,'ko-')
% hold on;
% plot(timeVals(endpoints2sample(SS_Bool)),SS_Sorted_Predict(SS_Bool),...
%     'ko','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% plot([time2SteadyStateInMinVal time2SteadyStateInMinVal], [0.5 1],'k--')
% xlabel('Time sampled (min)')
% ylabel('Prediction of sorting at steady state')
% set(gca,'XScale','log')
% 
% figure(600)
% plot(timeVals(endpoints2sample),Time_Sorted_Predict,'ko-')
% hold on;
% plot(timeVals(endpoints2sample(SS_Time_Bool)),Time_Sorted_Predict(SS_Time_Bool),...
%     'ko','MarkerFaceColor','k','LineWidth',1,'MarkerSize',3)
% plot([time2SteadyStateInMinVal time2SteadyStateInMinVal], 10.^[-2 8],'k--')
% xlabel('Time sampled (min)')
% ylabel('Time to sorting (min)')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% hold off;
% 
% % 
% % %%
% % 
% % PercentSortedWindow = 0.01;
% % time2SteadyStateInMin = interp1(percentSortedPredict,time2plotPredict,PercentSortedAtEq-PercentSortedWindow)
% % 
% % time2SteadyStateInHrs = time2SteadyStateInMin/60
% % 
% % time2SteadyStateInDays = time2SteadyStateInMin/(12*60)
% % 
% % time2SteadyStateInWeeks = time2SteadyStateInMin/(7*12*60)
% % 
% % time2SteadyStateInYears = time2SteadyStateInMin/(365*12*60)
% % 
% % %%
% % 
% % % 
% % % figure; plot(timeStepInMin_OverTime)