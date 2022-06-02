clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');


%%

color_gray = [0.3 0.3 0.3];

V = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model2/trace/%s/w%d/Vconnect.dat','weak',0.5*10));
%%

t1 = 0.01; 
t2 = 5000;
t=0.01:0.01:20000;
T1=t1/0.01; 
T2= t2/0.01;  

figure
plot(t(T1:T2), V(T1:T2,1),'color', color_gray, 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
xlim([108 2900])
ylim([-85 55])
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4 2]);
%print('/Users/kathleen/Documents/PhD/2022-NMOD/Fig_reset/zoomFig1_VE','-depsc')%, '-painters')



figure
plot(t(T1:T2), V(T1:T2,3),'color', 'k', 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', '-')
xlim([108 2900])
axis off
ylim([-85 55])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4 2]);
%print('/Users/kathleen/Documents/PhD/2022-NMOD/Fig_reset/zoomFig1_VC','-depsc')%, '-painters')


figure
plot(t(T1:T2), V(T1:T2,2),'color',color_gray, 'linewidth', 1)
xlim([108 2900])
ylim([-85 55])
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4 2]);
%print('/Users/kathleen/Documents/PhD/2022-NMOD/Fig_reset/zoomFig1_VI','-depsc')%, '-painters')

%%

t1 = 350; 
t2 = 850;
t=0.01:0.01:20000;
T1=t1/0.01; 
T2= t2/0.01;  

figure
plot(t(T1:T2), V(T1:T2,1),'color', color_gray, 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
xlim([t1 t2])
ylim([-85 55])
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3 1.5]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig2/zoomFig2_VE','-depsc')%, '-painters')



figure
plot(t(T1:T2), V(T1:T2,3),'color', color_gray, 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', '-')
xlim([t1 t2])
axis off
ylim([-85 55])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3 1.5]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig2/zoomFig2_VC','-depsc')%, '-painters')


%%

t1 = 1300; 
t2 = 3200;
t=0.01:0.01:20000;
T1=t1/0.01; 
T2= t2/0.01;  

figure
plot(t(T1:T2), V(T1:T2,1),'color', color_gray, 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
xlim([t1 t2])
ylim([-85 55])
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.7 1.5]);
print('/Users/kathleen/Documents/PhD/2022-NMOD/Fig_reset/zoomFig3_VE','-depsc')%, '-painters')



figure
plot(t(T1:T2), V(T1:T2,3),'color', color_gray, 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', '-')
xlim([t1 t2])
axis off
ylim([-85 55])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 3.7 1.5]);
print('/Users/kathleen/Documents/PhD/2022-NMOD/Fig_reset/zoomFig3_VC','-depsc')%, '-painters')

