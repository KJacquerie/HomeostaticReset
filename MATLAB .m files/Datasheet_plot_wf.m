clear all
close all
clc
pt=11;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%%
dt = 0.01;
Ttrans = 1000;
T = 50000;
t= dt:dt:T;


expm = 'F75';
type='lin'; 

f=[0.1; 1; 5; 10; 15; 20; 25; 30; 35; 40; 45;  50];

load('Sjo.mat')
%% Model 1  | Gamma

Gamma_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model1/results_tonic/checkF75_Gamma_wpos.dat');
Gamma_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model1/results_tonic/checkF75_Gamma_wneg.dat');
function_fig_wf(f, Gamma_wneg, Gamma_wpos, 'model1', 1.6)

%% Model  2 | Control

Control_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkF75_Control_wpos.dat');
Control_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkF75_Control_wneg.dat');
function_fig_wf(f, Control_wneg, Control_wpos, 'model2', 1.6)
%% Model 3 | DP (GB2106)
GB2016_DP_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model3/results_tonic/checkF75_kernel_DP_wpos.dat');
GB2016_DP_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model3/results_tonic/checkF75_kernel_DP_wneg.dat');
function_fig_wf(f, GB2016_DP_wneg, GB2016_DP_wpos, 'model3', 1.6)

%% Model 4
GB2012_cortex_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model4/results_tonic/checkF75_cortex_wpos.dat');
GB2012_cortex_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model4/results_tonic/checkF75_cortex_wneg.dat');
function_fig_wf(f, GB2012_cortex_wneg, GB2012_cortex_wpos, 'model4', 1.6)


%% Model 5
GB2012_DP_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model5/results_tonic/checkF75_DP_wpos.dat');
GB2012_DP_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model5/results_tonic/checkF75_DP_wneg.dat');
function_fig_wf(f, GB2012_DP_wneg, GB2012_DP_wpos, 'model5', 1.6)

%% Model 6 | Shouval
Shouval_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model6/results_tonic/checkF75_CC_wpos.dat');
Shouval_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model6/results_tonic/checkF75_CC_wneg.dat');
function_fig_wf(f, Shouval_wneg, Shouval_wpos, 'model6', 1.6)


%% Model 7 | Deperrois (dep)
Dep_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model7/results_tonic/checkF75_Dep_SJO_wpos.dat');
Dep_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model7/results_tonic/checkF75_Dep_SJO_wneg.dat');
function_fig_wf(f, Dep_wneg, Dep_wpos, 'model7',2)

%% Model 7 écrit en Omega
% Dep_DD_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model7/results_tonic/checkF75_Dep_DD_SJO_wpos.dat');
% Dep_DD_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model7/results_tonic/checkF75_Dep_DD_SJO_wneg.dat');
% function_fig_wf(f, Dep_DD_wneg, Dep_DD_wpos, 'model7_Omega',2)

%% Model 8 | Deperrois (STD)
STD_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_DD_SJO_wpos.dat');
STD_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_DD_SJO_wneg.dat');
function_fig_wf(f, STD_wneg, STD_wpos, 'model8', 2)

%% Model 8 | Deperrois (STD) écrit en Omega
STD_DD_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_DD_SJO_wpos.dat');
STD_DD_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_DD_SJO_wneg.dat');

function_fig_wf(f, STD_DD_wneg, STD_DD_wpos, 'model8_Omega', 2)


%% Model 9 (too fast)
%Shouval_DD_rel_wpos = load('/Users/kathleen/Documents/PhD/2022-NMOD/GBxShouval/results/checkF75_DD_rel_wpos.dat');
%Shouval_DD_rel_wneg = load('/Users/kathleen/Documents/PhD/2022-NMOD/GBxShouval/results/checkF75_DD_rel_wneg.dat');

%function_fig_wf(f, Shouval_DD_rel_wneg, Shouval_DD_rel_wpos, 'Shouval_DD_rel', 2)
%% Model 9 | Shouval DD rel 
Shouval_DD_rel_half_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model9/results_tonic/checkF60_DD_rel_half_wpos.dat');
Shouval_DD_rel_half_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model9/results_tonic/checkF60_DD_rel_half_wneg.dat');

function_fig_wf(f, Shouval_DD_rel_half_wneg, Shouval_DD_rel_half_wpos, 'model9', 2)


%% Model 10 | Shouval-CC-rel

Shouval_CC_rel_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model10/results_tonic/checkF60_CC_rel_wpos.dat');
Shouval_CC_rel_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model10/results_tonic/checkF60_CC_rel_wneg.dat');
function_fig_wf(f, Shouval_CC_rel_wneg, Shouval_CC_rel_wpos, 'model10', 2)


%% Model 11 = Model 7 en HB | Deperrois DD rel
Dep_DD_rel_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model11/results_tonic/checkFSJO_DD_rel_wpos.dat');
Dep_DD_rel_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model11/results_tonic/checkFSJO_DD_rel_wneg.dat');
function_fig_wf(f, Dep_DD_rel_wneg, Dep_DD_rel_wpos, 'model11',2)

%% Model 12
Dep_CC_rel_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model12/results_tonic/checkFSJO_CC_rel_wpos.dat');
Dep_CC_rel_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model12/results_tonic/checkFSJO_CC_rel_wneg.dat');
function_fig_wf(f, Dep_CC_rel_wneg, Dep_CC_rel_wpos, 'model12',2)


%% Model 13
STD_DD_rel_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model13/results_tonic/checkF60_STD_DD_rel_SJO_wpos.dat');
STD_DD_rel_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model13/results_tonic/checkF60_STD_DD_rel_SJO_wneg.dat');

function_fig_wf(f, STD_DD_rel_wneg, STD_DD_rel_wpos, 'model13',2)

%% Model 14
STD_CC_rel_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model14/results_tonic/checkF60_CC_rel_wpos.dat');
STD_CC_rel_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model14/results_tonic/checkF60_CC_rel_wneg.dat');
function_fig_wf(f, STD_CC_rel_wneg, STD_CC_rel_wpos, 'model14',2)




%%

figure
hold on
plot([0 50], [1 1 ], '-', 'color', [0.9 0.9 0.9], 'linewidth', 0.5)
plot(f,Control_wneg,':', 'color', 'k', 'linewidth', 1) 
plot(f,Control_wpos,'-','color', 'k', 'linewidth', 1) 
errorbar(Sjo(1:2:end,1),Sjo(1:2:end,3)+1, Sjo(1:2:end,4),'s','color', [0.5 0.5 0.5],'markersize',2.8)
errorbar(Sjo(2:2:end,1),Sjo(2:2:end,3)+1, Sjo(2:2:end,4),'o','color', [0.5 0.5 0.5],'markersize',2.8)
%le= legend('75','','SJO','');

xlim([-3 50])
ylim([0.25 2])
%set(le, 'fontsize',pt, 'interpreter','latex', 'location', 'northwest')
xlabel('$f [Hz]$', 'fontsize', pt, 'interpreter','latex')
ylabel('$\Delta w [-]$', 'fontsize', pt, 'interpreter','latex')
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7 4.5]);
%print('/Users/kathleen/Documents/PhD/2022-NMOD/Fig_wf/wf_Control_Cosyne','-depsc')

 
%%


figure
hold on
plot(f,Control_wneg,'-.', 'color', 'k', 'linewidth', 2.5) 
plot(f,Control_wpos,'-','color', 'k', 'linewidth', 2.5) 
errorbar(Sjo(1:2:end,1),Sjo(1:2:end,3)+1, Sjo(1:2:end,4),'s','color', [0.5 0.5 0.5],'markersize',10)
errorbar(Sjo(2:2:end,1),Sjo(2:2:end,3)+1, Sjo(2:2:end,4),'o','color', [0.5 0.5 0.5],'markersize',10)
%le= legend('75','','SJO','');
xlim([-3 50])
ylim([0.25 2])
%set(le, 'fontsize',pt, 'interpreter','latex', 'location', 'northwest')
xlabel('$f [Hz]$', 'fontsize', pt, 'interpreter','latex')
ylabel('$\Delta w [-]$', 'fontsize', pt, 'interpreter','latex')
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 13]);
%print('/Users/kathleen/Documents/PhD/2022-NMOD/Fig_wf/wf_GB2016_CosynePoster','-depsc')


