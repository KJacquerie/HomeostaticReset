close all 
clc
clear all 

dt1 = linspace(-80,80, 81);
dt50 = linspace(-10, 10, 21); 
pt=11;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

NB_models = 2; 
models=zeros(length(dt1),NB_models); 

%% Modele 3

models(:,1) = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model3/results_tonic/STDP_%s_f1.dat','kernel_DP'));

%% Modele 5
models(:,2) = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model5/results_tonic/STDP_%s_f1.dat','DP'));

%%

mean_models =mean(models,2)';
std_models = std(models,0,2)';
figure
hold on 
plot(dt1, models(:,1), 'color', [0.5 0.5 0.5],'linewidth', 0.5) 
plot(dt1, models(:,2), 'color', [0.5 0.5 0.5],'linewidth', 0.5) 
plot(dt1, mean_models, 'color', 'b','linewidth', 3)


figure
hold on 
errorbar(dt1, mean_models,std_models)
plot(dt1, mean_models, 'color', 'b','linewidth', 3)

curve1 = (mean_models + std_models);
curve2 = (mean_models - std_models);


figure
hold on
box off
plot([-100 100], [1 1], '-', 'color', [0.8 0.8 0.8], 'linewidth',0.5)
plot([0 0], [0 2], '-', 'color', [0.8 0.8 0.8], 'linewidth',0.5)
fill([dt1 fliplr(dt1)], [curve1  fliplr(curve2)], [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(dt1,mean_models, '-','color', 'k', 'linewidth', 1.25)

xlabel('$\Delta [ms]$', 'fontsize', pt, 'interpreter','latex')
ylabel('$\Delta w [-]$', 'fontsize', pt, 'interpreter','latex')

xlim([-80 80])
ylim([0.25 2])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 7 4.5]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig2/SUPERPOSITION_STDP','-depsc', '-painters')

print('/Users/kathleen/Documents/PhD/2022-reset/Fig2/SUPERPOSITION_STDP','-dsvg')%, '-painters')




