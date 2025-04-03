
%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This script sets the default simulation parameters

%% Set up the system

    % Clear the system
        close all; 
        clear all;

    % Choose where to store the parameters
        saveFilePath = ['\Garner_etal_2025_Fluidity_Sorting\Code\'...
            'Simulation_Code\basicModelParameters.mat'];

%% Choose the physical parameters

    % Choose the energies (in # kT_lab)
        E_het = 0;
        E_homo = -10^8;
        kT_eff = -E_homo;
    % Choose the size of the cell in nanometers
        cellRadius = 7*10^3;
    % Choose the distance between cell centers
        interNulcearDistance = 2*cellRadius;    
    % Choose the properties governing thermal energy in (picoNewton nanometers)
        % Boltzman constant in pN nm / K
            k_B = 0.0138;
        % Temperature in C
            T_C = 37;  
    % Choose the viscosity of the tissue in pN ms nm^{-2} (= 10^3 Pa s)
        etaVal = 1*10^(-3);

        etaValWater = 2.414*10^(-8)*10^(247.8/((273.15 + T_C)-140));
        etaValRel2Water = etaVal/etaValWater;

%% Choose the size of the simulation

    % Choose the number of cells along each dimension
        sizeX = 50;
        sizeY = 50;

    % Calculate the number of timepoints to save
        numTP2Save = 1000;
        % Check this is not over the max
            % Choose the maximum amount of data to store
                maxNumBytes = 10^9;
            % Choose the maximum number of timepoints to store
                maxNumTPs2Save = round(maxNumBytes./(sizeX*sizeY*8*3));
            if numTP2Save>maxNumTPs2Save
                numTP2Save = maxNumTPs2Save;
            end

%% Choose the simulation parameters

    % Choose the probability of an individual neighbor pair swapping
    % at each timepoint (lower = more accurate simulations)
        pSwap_0 = 0.005;

    % Choose the connectivity of the cells
    % (4 for 4- or 8 for 8-connected neighborhoods)
        connectivity = 4;

    % Choose the boundary conditions (0 for periodic or 1 for fixed)
        BCs = 1;
     
%% Save the file

    % Save all parameters to a structure and clear all of the remaining variables
        inputVars = v2struct;
        % Clear all other variables
            clearvars -except inputVars

    % Save the file
        save(inputVars.saveFilePath,'-struct','inputVars','-v7.3')
