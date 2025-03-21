function [x] = evalInverseAsymHill(a,b,c,d,y)

%% The following code was written by Rikki Garner to generate the figures in 
% Tissue Fluidity: A Double-Edged Sword for Multicellular Patterning
% Rikki M. Garner, Sean E. McGeary, Allon M. Klein, Sean G. Megason
% bioRxiv 2025.03.01.640992; doi: https://doi.org/10.1101/2025.03.01.640992
% This code was last updated on 2025/3/20

% This script fits data to an inverse asymmetric hill function


x = b*(((a-0.5)/(y-0.5))^(1/d)-1)^(-1/c);