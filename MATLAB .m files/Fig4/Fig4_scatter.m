clear all
close all
clc

pt=11;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%
scatter_simuvspredict('Wake', 'small1');
%%
scatter_simuvspredict('Sleep', 'small2');
%%
scatter_simuvspredict('Wake', 'large');%%
