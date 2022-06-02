close all 
clc
clear all 

dt1 = linspace(-80,80, 81);
dt50 = linspace(-10, 10, 21); 
pt=11;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


%%  Model 2

STDP_Control_f1 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/STDP_%s_f1.dat','Control'));
STDP_Control_f50 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/STDP_%s_f50.dat','Control'));


function_fig_STDP(dt1,1,STDP_Control_f1,'model2')
function_fig_STDP(dt50,50,STDP_Control_f50,'model2')
%% Modele 3

STDP_kernel_DP_f1 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model3/results_tonic/STDP_%s_f1.dat','kernel_DP'));
STDP_kernel_DP_f50 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model3/results_tonic/STDP_%s_f50.dat','kernel_DP'));


function_fig_STDP(dt1,1,STDP_kernel_DP_f1,'model3')
function_fig_STDP(dt50,50,STDP_kernel_DP_f50,'model3')

%% Model 4
STDP_GB2012_cortex_f1 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model4/results_tonic/STDP_%s_f1.dat','cortex'));
STDP_GB2012_cortex_f50 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model4/results_tonic/STDP_%s_f50.dat','cortex'));


function_fig_STDP(dt1,1,STDP_GB2012_cortex_f1,'model4')
function_fig_STDP(dt50,50,STDP_GB2012_cortex_f50,'model4')

%% Modele 5
STDP_GB2012_DP_f1 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model5/results_tonic/STDP_%s_f1.dat','DP'));
STDP_GB2012_DP_f50 = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model5/results_tonic/STDP_%s_f50.dat','DP'));


function_fig_STDP(dt1,1,STDP_GB2012_DP_f1,'model5')
function_fig_STDP(dt50,50,STDP_GB2012_DP_f50,'model5')

