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

%% Model  2 | Control - 75 pulses

Control_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkF75_Control_wpos.dat');
Control_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkF75_Control_wneg.dat');


%% Model  2 | Control - 5x15 pulses

SJO_wpos = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkFSJO_Control_wpos.dat');
SJO_wneg = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_tonic/checkFSJO_Control_wneg.dat');

%%
pt=11;
figure
hold on
plot([0 50], [1 1 ], '-', 'color', [0.9 0.9 0.9], 'linewidth', 0.5)
plot(f,Control_wneg,':', 'color', 'k', 'linewidth', 1) 
plot(f,Control_wpos,'-','color', 'k', 'linewidth', 1) 
plot(f,SJO_wneg,':', 'color', [0. 0.4 0.4], 'linewidth', 1) 
plot(f,SJO_wpos,'-','color', [0. 0.4 0.4], 'linewidth', 1) 
errorbar(Sjo(1:2:end,1),Sjo(1:2:end,3)+1, Sjo(1:2:end,4),'s','color', [0.5 0.5 0.5],'markersize',2.8)
errorbar(Sjo(2:2:end,1),Sjo(2:2:end,3)+1, Sjo(2:2:end,4),'o','color', [0.5 0.5 0.5],'markersize',2.8)
%le= legend('75','','SJO','');

xlim([-1 51])
ylim([0.4 2])
%ylim([0.4 2])
%set(le, 'fontsize',pt, 'interpreter','latex', 'location', 'northwest')
xlabel('$f [Hz]$', 'fontsize', pt, 'interpreter','latex')
ylabel('$\Delta w [-]$', 'fontsize', pt, 'interpreter','latex')

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7 4.5]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig_SI/75vsSJO', '-depsc');

