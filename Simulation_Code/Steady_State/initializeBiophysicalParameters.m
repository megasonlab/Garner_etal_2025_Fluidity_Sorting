function [inputVars] = initializeBiophysicalParameters(inputVars)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by setUpAndRun_Sorting*.m, in order to initialize
% the biophysical parameters

    % Unpack this structure and clear the original structure
        v2struct(inputVars)
        clear inputVars        
            
    % Calculate the diffusion-limitedtransition rate
        % Calculate the thermal energy in (picoNewton nanometers)
            % Calcuate the temperature in K
                T_K = 273.15 + T_C;
            % kT   
                kT = k_B*T_K;
        % Calculate the Stokes drag for this segment size (assuming the dynamic
        % viscoty of water (in picoNewtons milliseconds / nanometer)
            gammaVal = 6*pi*etaVal*cellRadius;
        % Calculate the rate of neighbor exchange per neighbor per second
            % 4*4  = your target could swap with any of their 4 neighbors, 
            % and each of your 4 neighbors could switch with any of their outter neighbors, 
            rSwap_0 = kT_eff*kT*10^3/(gammaVal*interNulcearDistance^2);
    
    % Load these variables into a structure
        inputVars = v2struct();

end
