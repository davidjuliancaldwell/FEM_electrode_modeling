% This is a script to viusalize the results of the FEM modeling
%
% David J. Caldwell 
%
% Requires: 
% Flux Matrix - m x n x p matrix
% Volt matrix - m x n x p matrix 
%
%% clear workspace, load data 
clear all;clc

% load in the different layers
load('FuxMagTest1.mat')
load('VoltTest1.mat')
%%

% look at Volt
threshold = 0;
plot_results_cross_section(VoltMatrix,threshold)

% look at Flux 
threshold = 0;
plot_results_cross_section(FluxMatrix,threshold)

%%
PlotSlices_iterate_results