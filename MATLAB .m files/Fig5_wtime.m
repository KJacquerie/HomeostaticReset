clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');




%%

function_plot_wtime_NMOD('calcium','Ampl','cortex','calcium_ampl',  3) % #9
function_plot_wtime_NMOD('calcium','tau','cortex','calcium_tau',  3) % #9

function_plot_wtime_TAG('calcium','Ampl_TAG','cortex','calcium_ampl_tag',  3) % #9
function_plot_wtime_TAG('calcium','tau_TAG','cortex','calcium_tau_tag',  3) % #9


%%
function_plot_wtime_NMOD('phenom','Ampl','Pair_SB','phenom_ampl',  3) % #9
function_plot_wtime_NMOD('phenom','tau','Pair_SB','phenom_tau',  3) % #9

function_plot_wtime_TAG('phenom','Ampl_TAG','Pair_SB','phenom_ampl_tag',  3) % #9
function_plot_wtime_TAG('phenom','tau_TAG','Pair_SB','phenom_tau_tag',  3) % #9


%%%%%%%%%

%%
close all
function_plot_wtimeFENS('calcium/normal','cortex','calcium_normal',  3) % #9
%%
close all
function_plot_wtimeFENS('calcium/omega','cortex','calcium_omega',3) % #1
%%
close all
function_plot_wtimeFENS('calcium/Ampl','cortex','calcium_ampl',  3) % #9

%%
close all
function_plot_wtimeFENS('calcium/Ampl_TAG','cortex','calcium_ampl_tag',3) % #1

%%
close all
function_plot_wtimeFENS('calcium/tau','cortex','calcium_tau',  3) % #9

%%
close all
function_plot_wtimeFENS('calcium/tau_TAG','cortex','calcium_tau_tag',  3) % #9

%%
close all
function_plot_wtimeFENS('calcium/tau_TAG','hippocampus','calcium_tau_hipp',  3) % #9


%%
close all
function_plot_wtimeFENS('phenom/normal','Pair_SB','Pair_normal',6) % #1

%%
close all
function_plot_wtimeFENS('phenom/Ampl','Pair_SB','Pair_ampl',6) % #1

%%
close all
function_plot_wtimeFENS('phenom/Ampl_TAG','Pair_SB','Pair_ampl_tag',6) % #1

%%
close all
function_plot_wtimeFENS('phenom/tau','Pair_SB','Pair_tau',6) % #1

%%
close all
function_plot_wtimeFENS('phenom/tau_TAG','Pair_SB','Pair_tau_tag',6) % #1



%%
% close all
% function_plot_wtime('omega','cortex','omega_cortex',6) % #1
% %%
% close all
% function_plot_wtime('omega','hippocampus','omega_hipp',  3) % #9
% %%
% close all
% function_plot_wtime('omega_winit','cortex','omegawinit_cortex',6) % #1
% %%
% close all
% function_plot_wtime('omega','hippocampus','omegawinit_hipp',  3) % #9
% 
% %%
% close all
% function_plot_wtime('tau','hippocampus','tau_hipp',  3) % #9