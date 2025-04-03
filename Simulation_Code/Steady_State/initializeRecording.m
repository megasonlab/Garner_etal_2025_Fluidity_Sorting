function [inputVars] = initializeRecording(inputVars)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by setUpAndRun_Sorting*.m, in order to initialize
% the temporal recording parameters

    % Unpack this structure and clear the original structure
        v2struct(inputVars)
        clear inputVars        
            
    % Initialize the time and counter variables
        timeSimInMin = 0;
        tp = 0;

    % Chose the timepoints to save the data
        lastTimeStep2Sample = 10000;
        timesSteps2Sample = unique(round(10.^(linspace(0,log10(lastTimeStep2Sample),numTP2Save-1))));
        nextTimeStep2Sample = timesSteps2Sample(1);
        tp2Record = 2;

    % Choose the timepoint to check if simulation has reached steady state
        nextTimeInMin2CheckSS = (pSwap_0/rSwap_0/60)*1000;

    % Pre-allocate the steady state boolean
        notSteadyState = true;

    % Load these variables into a structure
        inputVars = v2struct();

end
