clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%
close all
function_plot_wtime_Article('model2','sleep','_RESET', 'weak', 6) % #1

close all
function_plot_wtime_Article('model9','sleep','_SATURATION', 'weak', 3) % #9
