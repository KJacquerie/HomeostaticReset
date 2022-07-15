clear all
close all
clc

pt=11;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

nPre = 100; 
nPost = 100; 

color_blue = [106 153 208]./255; 

%%

display_correlation('Wake', 'large');
% display_correlation('Sleep',  'large');

