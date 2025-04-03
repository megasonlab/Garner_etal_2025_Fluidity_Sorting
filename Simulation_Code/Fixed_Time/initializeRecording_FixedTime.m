function [inputVars] = initializeRecording_FixedTime(inputVars)


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

    % Choose the timepoints to sample
        % % Linear
           % times2SampleInMin = 0:(time2SimInMin/(numTP2Save-1)):time2SimInMin;
        % Log-distributed
            firstTime2Sample = 10^(-4);
            lastTime2Sample = time2SimInMin;
            interval = (log10(lastTime2Sample)-log10(firstTime2Sample))/(numTP2Save-2);
            times2SampleInMin = [0 10.^(log10(firstTime2Sample):interval:log10(lastTime2Sample))];

    % Initialize the time and counter variables
        timeSimInMin = 0;
        tp = 0;

    % Chose the timepoints to save the data
        nextTime2SampleInMin = times2SampleInMin(2);
        tp2Record = 2;                

    % Load these variables into a structure
        inputVars = v2struct();

end
