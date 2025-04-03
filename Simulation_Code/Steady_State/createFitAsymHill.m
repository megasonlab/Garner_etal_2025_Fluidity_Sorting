function [fitresult, gof] = createFitAsymHill(xVals, yVals, plotResults)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by simSorting_FeedFile_SS.m, and fits the sorting
% curve to an asymmetric hill function in order to
% determine whether the simulation has reached steady state


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xVals, yVals );

% Set up fittype and options.
ft = fittype( '0.5+(a-0.5)/(1+((b/x)^c))^d', 'independent', 'x', 'dependent', 'y' );
excludedPoints = (xData <= NaN);%&(xData <= 10^(-2));
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.5 0 0 0.3];
opts.StartPoint = [0.25 0.225921780972399 0.5 1];
opts.Upper = [1 Inf Inf Inf];
opts.Exclude = excludedPoints;
%opts.Weights = exp(-(1:length(yData))/100);

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if plotResults
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData, excludedPoints );
legend( h, 'yVals vs. xVals', 'Excluded yVals vs. xVals', 'untitled fit 1', 'Location', 'NorthWest', 'Interpreter', 'none' );
% Label axes
xlabel( 'xVals', 'Interpreter', 'none' );
ylabel( 'yVals', 'Interpreter', 'none' );
set(gca,'XScale','log')
grid on
drawnow;
end


