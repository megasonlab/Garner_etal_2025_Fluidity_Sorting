function [x] = evalInverseAsymHill(a,b,c,d,y)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/4/3

% This function is called by simSorting_FeedFile_SS.m, and evaluates the 
% inverse asymmetric hill function in order to
% determine whether the simulation has reached steady state

x = b*(((a-0.5)/(y-0.5))^(1/d)-1)^(-1/c);