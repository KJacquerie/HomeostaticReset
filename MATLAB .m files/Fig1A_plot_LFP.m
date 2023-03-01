clear all
close all
clc 



color_blue = [106 153 208]./255; 
color_gray = [ 0.3 0.3 0.3]; 
color_pink = [250/255 244/255 247/255];

color_corr = [63/255 92/255 206/255]; 
color_uncorr = [130/255 187/255 255/255]; 

%%
directory_name = "/Users/kathleen/Documents/PhD/2022-reset"; 
% time vector
T = 60000;
dt = 0.01;
Tdt = T/dt;
t = dt:dt:T;


%%

% LFP vector
LFP_E= load(sprintf('%s/fig1A_LFP/results_LFP/LFP_E.dat',directory_name));
LFP_I= load(sprintf('%s/fig1A_LFP/results_LFP/LFP_I.dat',directory_name));
LFP_C= load(sprintf('%s/fig1A_LFP/results_LFP/LFP_C.dat',directory_name));
%%
V= load(sprintf('%s/LFP/results_LFP/V.dat',directory_name));

%%
T1=13500;
T2=16900;
interv = T1/dt:T2/dt;

yy=smooth(LFP_C,2500);
figure
subplot(2,1,1)
plot(t(interv)*1e-3, LFP_C(interv)) 
subplot(2,1,2)
plot(t(interv)*1e-3, yy(interv)) 
%%
tshort = T1:T2;
figure
%subplot(4,1,1)
plot(tshort,  V(6,T1:T2)', 'color', color_gray, 'linewidth', 0.75)
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
xlim([T1-100 T2+100])
box off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5.7 0.95]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig_LFP/Fig1A_Vpre','-depsc')%, '-painters')

%%
figure
%subplot(4,1,2)
plot(tshort,  V(12,T1:T2)', 'color', color_gray, 'linewidth', 0.75)
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
xlim([T1 T2])
box off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5.7 0.95]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig_LFP/Fig1A_Vpost','-depsc')%, '-painters')

%%
figure
%subplot(4,1,3)
plot(t(interv)*1e-3, LFP_C(interv)) 

%%
figure
%subplot(4,1,4)
plot(t(interv)*1e-3, yy(interv), 'color', color_gray, 'linewidth', 0.75) 
ylim([-0.05 0.76])
xlim([(T1-100)/1000 (T2+100)/1000])
yticks([0 0.5])
yticklabels({'',''})
xticks([0])
xticklabels({''})

box off
%axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5.7 1.05]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig_LFP/Fig1A_LFP','-depsc')%, '-painters')


%%
noverlap = 100/dt; %100
window = 10000/dt; %5000
fs=1000/dt;

spectrogram(yy,window,noverlap,[],fs,'MinThreshold', -105,'yaxis')
ylim([0 40*1e-3])


%%
% %% ZOOM FIGURE
% interv = 20500/dt:23000/dt; 
% pt=11;
% ptx = 6.2; 
% pty = 1.75; 
% color_vec = [1/256 1/256 1/256]*110;
% 
% 
% figure
% plot(t(interv)*1e-3, LFP_I(interv), 'color', color_vec)
% xlim([min(t(interv))*1e-3 max(t(interv))*1e-3])
% box off
% title('')
% xticks([21,22])
% xticklabels({'0','1'})
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',pt)
% set(gcf,'PaperPositionMode','auto');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 ptx pty]);
% % print(sprintf('Figures/V/LFPtrace%s_%d_zoom', Name, var),'-depsc')
% % print(sprintf('Figures/V/LFPtrace%s_%d_zoom', Name, var),'-dpdf')
% 
% %%