clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%
close all
function_plot_wtime_Oneline('model2','sleep','_RESET', 'weak', 6) % #1
%%
close all
function_plot_wtime_Oneline('model9','sleep','_SATURATION', 'weak', 1) % #9
