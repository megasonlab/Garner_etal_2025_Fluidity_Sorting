function [meanK_2DAcorr_2DFFT] = calculateDomainSizeFromStructureFactor(cellType_onGrid_OverTime,timepoints2Analyze,sizeX,doPlot)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This function takes in a 2D cell type matrix and calulcates the 
% mean wavenumber from the structure factor over time

% Pre-allocate space to store the mean wavenumber
    meanK_2DAcorr_2DFFT = nan([length(timepoints2Analyze) 1]);

% Initialize the counter
    tpCounter = 0;

% Loop through each timepoint and calculate the mean domain size
for numTP = timepoints2Analyze

    % Update the counter
        tpCounter = tpCounter+1;

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
            meanK_2DAcorr_2DFFT(tpCounter) = sum(waveNumbers.*cellType_2Dacorr_2DFFT_RMean')./sum(cellType_2Dacorr_2DFFT_RMean);

    % Plot the steps
    if doPlot
        figure(1991)
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
        plot([meanK_2DAcorr_2DFFT(tpCounter) meanK_2DAcorr_2DFFT(tpCounter)],yl,'m--','LineWidth',2)
        hold off;
        legend({'',sprintf('Mean k = %2.2f',meanK_2DAcorr_2DFFT(tpCounter))})
        title({'Radially-averaged 2D FFT of','2D spatial autocorrelation function','S(k,t) = FFT(C(r,t))','(equal time structure factor)'})
        ylabel('Mean FFT')
        xlabel('Wavenumber (k)')
        drawnow
        drawnow
    end
end

end