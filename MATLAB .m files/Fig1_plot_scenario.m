clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

%%
% idx from 1 to 3 corresponds to correlated inputs > idx=2
% idx from 4 to 6 corresponds to uncorrelated inputs > idx = 6
% to see different plots change the idx
idx=6;


for idx_w=1:1:6
    w(:,idx_w) = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/scenario/w_neuron%d.dat', idx_w));
end


color_blue = [106 153 208]./255; 
color_gray = [ 0.3 0.3 0.3]; 
color_pink = [250/255 244/255 247/255];

color_corr = [63/255 92/255 206/255]; 
color_uncorr = [130/255 187/255 255/255]; 
%%

t1 = 0.01; 
t2 = 160000;
t=0.01:1:160000;
T1=t1/0.01; 
T2= t2/0.01;  

figure
hold on 
v = [40000 0; 80000 0;80000 1; 40000 1];
v2 = [120000 0; 160000 0;160000 1; 120000 1];
f = [1 2 3 4];
patch('Faces',f,'Vertices', v, 'FaceColor', color_pink, 'EdgeColor', 'none')

patch('Faces',f,'Vertices', v2, 'FaceColor', color_pink, 'EdgeColor', 'none')
plot(t, w(:,1:3), 'color', color_corr, 'linewidth', 1)
plot(t, w(:,4:6), 'color', color_uncorr, 'linewidth', 1)
box off
ylim([0 1])
%xticks([0 1000])
%xticklabels({'0','1000'})
%yticks([0 1])
%yticklabels({'0','1'})
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 10 4]);
print('/Users/kathleen/Documents/PhD/2022-reset/Fig_scenario/trace_scenario','-depsc')%, '-painters')



%% plot switch figure 1A

tstart = 38600;
tstop=41600;

expm = 'uncorr'; 
VE = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/scenario/V_E%d.dat', 5));
VC = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/scenario/V_C%d.dat', 5));

figure
plot(t(tstart:tstop), VE(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-100 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig_scenario/Fig1A_VE_%s',expm),'-depsc')%, '-painters')



figure
plot(t(tstart:tstop), VC(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-100 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig_scenario/Fig1A_VC_%s',expm),'-depsc')%, '-painters')




%%

VE = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/scenario/V_E%d.dat', idx));
VC = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/scenario/V_C%d.dat', idx));

%%
if(idx<=3) 
    expm='corr'; 
    color_plot = color_corr; 
else
    expm='uncorr'; 
    color_plot = color_uncorr;
end

%% plot the switch associated to the desired correlated curve
tstart = 38600;
tstop=41600;

figure
plot(t(tstart:tstop), VE(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-100 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/VE_%s',expm),'-depsc')%, '-painters')
%%

figure
plot(t(tstart:tstop), VC(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-100 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/VC_%s',expm),'-depsc')%, '-painters')




%% plot the zoom tonic
tstart = 100;
tstop=1100;
zoom = 'T';

figure
plot(t(tstart:tstop), VE(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-50 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/VE_%s%s',expm,zoom),'-depsc')%, '-painters')


figure
plot(t(tstart:tstop), VC(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-50 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/VC_%s%s',expm, zoom),'-depsc')%, '-painters')
%%
figure
plot(t(tstart:tstop), w(tstart:tstop, idx), 'color',color_plot, 'linewidth', 1)
xlim([tstart-50 tstop])
box off
xticks(0)
xticklabels({''})
yticks([ 0])
yticklabels({''})
if(idx<=3)
    ylim([0.499 0.51])
else
    ylim([0.4895 0.5005])
end
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/w_%s%s',expm, zoom),'-depsc')%, '-painters')


%% plot the zoom burst
tstart = 125510;
tstop=126510;
zoom = 'B';

figure
plot(t(tstart:tstop), VE(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-50 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/VE_%s%s',expm,zoom),'-depsc')%, '-painters')


figure
plot(t(tstart:tstop), VC(tstart:tstop), 'color',color_gray, 'linewidth', 1)
xlim([tstart-50 tstop])
box off
xticks(0)
xticklabels({''})
yticks([-50 0])
yticklabels({'',''})
ylim([-90 25])
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/VC_%s%s',expm, zoom),'-depsc')%, '-painters')
%%
figure
plot(t(tstart:tstop), w(tstart:tstop,idx), 'color',color_plot, 'linewidth', 1)
xlim([tstart-50 tstop])
box off
xticks(0)
xticklabels({''})
yticks([ 0])
yticklabels({''})
if(idx<=3)
    ylim([0.64 0.6602 ])
else
    ylim([0.44 0.46])
end
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 1]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig1_scenario/w_%s%s',expm, zoom),'-depsc')%, '-painters')



