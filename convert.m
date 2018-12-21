%% Save output and convert to m and c code fragments
clear;
clc;
% load equations for predictions and updates
load('PredictionEquations.mat');

fileName = strcat('SymbolicOutput',int2str(nStates),'.mat');
save(fileName);
SaveScriptCode(nStates);
ConvertToM(nStates);
ConvertToC(nStates);