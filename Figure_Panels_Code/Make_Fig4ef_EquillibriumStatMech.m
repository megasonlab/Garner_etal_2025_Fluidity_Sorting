
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

%% Fig. 4e -- Plot the probability of being sorted as a function of the differential adhesion energy 

    % Choose the energy difference between the adhesion energy in
    % the sorted state (E_S) and the anti-sorted state (E_A)
        dE_Vals = [-2:0.01:2];
    % Choose the thermal energies
        kT_Vals = [0.1 0.3 1 3 100];  
    % Choose the colormap
        CM = jet(length(kT_Vals))*0.8;
        CM(4,:) = CM(5,:);
        CM(5,:) = [1 0 0]*0.8;

    % Pull up the figure
        figure(1)
    % Plot the sorting probabilty vs differential adhesion energy for
    % different choices of kT
        for kT_Num = 1:length(kT_Vals)    
        plot(dE_Vals,...
            1./(1+(0.5*exp((dE_Vals)./ ...
                kT_Vals(kT_Num)))),'--','LineWidth',2,...
                'Color',[CM(kT_Num,:) 0.5])
        hold on;
        end
        hold off;
    % Label the plot
        ylim([0 1])
        xlabel('Differential adhesion energy (E_S - E_U)')
        ylabel('Probability of being in the sorted state')
        title('Varying E_M')
        legend([cellstr(num2str(kT_Vals', 'E_M=%1.1f'))],...
            'Location','Southwest')
        set(gca,'FontName','Helvetica','FontSize',5); 

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig4e_ProbabilityOfSortingVsdEA'];
        % Choose the figure size and position
            figHandle = figure(1);
            figHandle.Position =  [475   325   210   160];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])

%% Fig. 4f -- Plot the probability of being sorted as a function of the motility energy (E_M)

    % Choose the energy difference between the adhesion energy in
    % the sorted state (E_S) and the anti-sorted state (E_U)
        dE_Vals = -2;
    % Choose the thermal energies
        kT_Vals = 10.^(linspace(-1,log10(20),100));  
  
    % Pull up the figure
        figure(1)
    % Plot the sorting probabilty vs differential adhesion energy for
    % different choices of kT  
        plot(kT_Vals,...
            1./(1+(0.5*exp((dE_Vals)./ ...
                kT_Vals))),'--','LineWidth',2,...
                'Color',[0 1 0]*0.9)
        hold on;
    % Label the plot
        ylim([0.65 1])
        xlim([kT_Vals(1) kT_Vals(end)])
        xl = xlim();
        plot(xl,0.67*[1 1],'k--')
        hold off;
        xlabel('E_M')
        ylabel('p_S')
        legend(sprintf('(E_S-E_U)=%1.1d',dE_Vals),'Entropic limit (random)')
        set(gca,'FontName','Helvetica','FontSize',5); 

     % Save the figure 
        % Choose the filename and path for the figure
            destinationTrack = [figOutputFolderPath 'Fig4f_ProbabilityOfSortingVsEM'];
        % Choose the figure size and position
            figHandle = figure(1);
            figHandle.Position =  [475   325   210   160];
        % Save the file as a pdf
            exportgraphics(gcf,[destinationTrack '.pdf'],'ContentType','vector')
        % Save as a png                
            saveas(gcf,[destinationTrack '.png'])
