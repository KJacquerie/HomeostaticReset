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

NB_models=14; 
m_wpos = zeros(length(f),NB_models); 
m_wneg = zeros(length(f),14);


%% Model 1

m_wpos(:,1) = load('/Users/kathleen/Documents/PhD/2022-reset/model1/results_tonic/checkF75_Gamma_wpos.dat');
m_wneg(:,1) = load('/Users/kathleen/Documents/PhD/2022-reset/model1/results_tonic/checkF75_Gamma_wneg.dat');

%% Model 2
m_wpos(:,2) = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkF75_Control_wneg.dat');
m_wneg(:,2) = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkF75_Control_wneg.dat');
 

%% Model 3
%m_wpos(:,3) = load('/Users/kathleen/Documents/PhD/2022-NMOD/GB_kernel/results/checkF75_kernel_DP_wpos.dat');
%m_wneg(:,3) = load('/Users/kathleen/Documents/PhD/2022-NMOD/GB_kernel/results/checkF75_kernel_DP_wneg.dat');


%% Model 4
m_wpos(:,3) = load('/Users/kathleen/Documents/PhD/2022-reset/model4/results_tonic/checkF75_cortex_wpos.dat');
m_wneg(:,3) = load('/Users/kathleen/Documents/PhD/2022-reset/model4/results_tonic/checkF75_cortex_wneg.dat');


%% Model 5
%m_wpos(:,5) = load('/Users/kathleen/Documents/PhD/2022-NMOD/GB2012/results/checkF75_DP_wpos.dat');
%m_wneg(:,5) = load('/Users/kathleen/Documents/PhD/2022-NMOD/GB2012/results/checkF75_DP_wneg.dat');

%% Model 6
m_wpos(:,4) = load('/Users/kathleen/Documents/PhD/2022-reset/model6/results_tonic/checkF75_CC_wpos.dat');
m_wneg(:,4) = load('/Users/kathleen/Documents/PhD/2022-reset/model6/results_tonic/checkF75_CC_wneg.dat');


%% Model 7

m_wpos(:,5) = load('/Users/kathleen/Documents/PhD/2022-reset/model7/results_tonic/checkF75_Dep_SJO_wpos.dat');
m_wneg(:,5) = load('/Users/kathleen/Documents/PhD/2022-reset/model7/results_tonic/checkF75_Dep_SJO_wpos.dat');

%% Model 8
%STD_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_SJO_wpos.dat');
%STD_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_SJO_wneg.dat');


%% Model 8 (écrit en forme standard) 
m_wpos(:,6) = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_DD_SJO_wpos.dat');
m_wneg(:,6) = load('/Users/kathleen/Documents/PhD/2022-reset/model8/results_tonic/checkF75_STD_DD_SJO_wneg.dat');


%% Model 9 
m_wpos(:,7) = load('/Users/kathleen/Documents/PhD/2022-reset/model9/results_tonic/checkF60_DD_rel_half_wpos.dat');
m_wneg(:,7) = load('/Users/kathleen/Documents/PhD/2022-reset/model9/results_tonic/checkF60_DD_rel_half_wneg.dat');



%% Model 10

m_wpos(:,8) = load('/Users/kathleen/Documents/PhD/2022-reset/model10/results_tonic/checkF60_CC_rel_wpos.dat');
m_wneg(:,8) = load('/Users/kathleen/Documents/PhD/2022-reset/model10/results_tonic/checkF60_CC_rel_wneg.dat');


%% Model 11 = Model 7 en HB
m_wpos(:,9) = load('/Users/kathleen/Documents/PhD/2022-reset/model11/results_tonic/checkFSJO_DD_rel_wpos.dat');
m_wneg(:,9) = load('/Users/kathleen/Documents/PhD/2022-reset/model11/results_tonic/checkFSJO_DD_rel_wneg.dat');


%% Model 12
m_wpos(:,10) = load('/Users/kathleen/Documents/PhD/2022-reset/model12/results_tonic/checkFSJO_CC_rel_wpos.dat');
m_wneg(:,10) = load('/Users/kathleen/Documents/PhD/2022-reset/model12/results_tonic/checkFSJO_CC_rel_wneg.dat');


%% Model 13
m_wpos(:,11) = load('/Users/kathleen/Documents/PhD/2022-reset/model13/results_tonic/checkF60_STD_DD_rel_SJO_wpos.dat');
m_wneg(:,11) = load('/Users/kathleen/Documents/PhD/2022-reset/model13/results_tonic/checkF60_STD_DD_rel_SJO_wneg.dat');


%% Model 14
m_wpos(:,12) = load('/Users/kathleen/Documents/PhD/2022-reset/model14/results_tonic/checkF60_CC_rel_wpos.dat');
m_wneg(:,12) = load('/Users/kathleen/Documents/PhD/2022-reset/model14/results_tonic/checkF60_CC_rel_wneg.dat');


%%

figure
hold on
plot([0 50], [1 1 ], '-', 'color', [0.9 0.9 0.9], 'linewidth', 0.5)
plot(f,m_wneg,':', 'color', 'b', 'linewidth', 1) 
plot(f,m_wpos,'-','color', 'r', 'linewidth', 1) 
plot(f, mean(m_wneg,2), 'color', 'b','linewidth', 3) 
plot(f, mean(m_wpos,2), 'color', 'r','linewidth', 3) 
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
%print('/Users/kathleen/Documents/PhD/2022-reset/Fig2/SUPERPOSITION','-depsc')

 
%%


figure
plot(f, mean(m_wneg,2), 'color', 'b','linewidth', 3) 
plot(f, mean(m_wpos,2), 'color', 'r','linewidth', 3) 

%%
mean_neg = mean(m_wneg,2)'; % your mean vector;
std_neg = std(m_wneg,0,2)';
curve1_neg = (mean_neg + std_neg);
curve2_neg = (mean_neg - std_neg);

mean_pos = mean(m_wpos,2)'; % your mean vector;
std_pos = std(m_wpos,0,2)';
curve1_pos = (mean_pos + std_pos);
curve2_pos = (mean_pos - std_pos);


figure
hold on 
errorbar(f, mean_neg, std_neg)
errorbar(f, mean_pos, std_pos)

%%
figure
hold on
box off
fill([f' fliplr(f')], [curve1_pos  fliplr(curve2_pos)], [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

fill([f' fliplr(f')], [curve2_neg  fliplr(curve1_neg)], [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(f,mean_neg, ':','color', 'k', 'linewidth', 1.25)
plot(f,mean_pos, '-','color', 'k', 'linewidth', 1.25)

errorbar(Sjo(1:2:end,1),Sjo(1:2:end,3)+1, Sjo(1:2:end,4),'s','color', [0.5 0.5 0.5],'markersize',2.8)
errorbar(Sjo(2:2:end,1),Sjo(2:2:end,3)+1, Sjo(2:2:end,4),'o','color', [0.5 0.5 0.5],'markersize',2.8)
xlabel('$f [Hz]$', 'fontsize', pt, 'interpreter','latex')
ylabel('$\Delta w [-]$', 'fontsize', pt, 'interpreter','latex')

xlim([-3 50])
ylim([0.25 2])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7 4.5]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig2/SUPERPOSITION_wf','-depsc', '-painters')


print('/Users/kathleen/Documents/PhD/2022-reset/Fig2/SUPERPOSITION_wf','-dsvg')%, '-painters')

