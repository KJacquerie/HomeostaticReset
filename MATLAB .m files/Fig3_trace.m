clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

ptx = 2.25;
pty = 0.75;
%% TRACE FIG 3A

for idx_Ivec = 1:1:8
V = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model2/trace/burst/Vconnect_%d.dat',idx_Ivec));


t1 = 0.01; 
t2 = 5000;
t=0.01:0.01:t2;
T1=4100/0.01; 
T2= 4900/0.01;  

figure
plot(t(T1:T2), V(T1:T2,1),'color', [0.3 0.3 0.3], 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
box off
ylim([-85 55])
%ylabel('$V_E$ [mV]', 'interpreter', 'latex', 'fontsize', 11)
yticks([-50 0 50])
yticklabels({'', '',''})
xticks([2000 3000 4000 5000])
xticklabels({'','','',''})
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig3_trace/Fig3_VE_I%d',idx_Ivec),'-depsc')%, '-painters')


figure
plot(t(T1:T2), V(T1:T2,3),'color', [0.3 0.3 0.3], 'linewidth', 1)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
box off
ylim([-85 55])
%ylabel('$V_E$ [mV]', 'interpreter', 'latex', 'fontsize', 11)
yticks([-50 0 50])
yticklabels({'', '',''})
xticks([2000 3000 4000 5000])
xticklabels({'','','',''})
axis off
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 ptx pty]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig3_trace/Fig3_VC_I%d',idx_Ivec),'-depsc')%, '-painters')
end


%%


VI2 = load('/Users/kathleen/Documents/PhD/2022-reset/model2/results_Variab/Vconnect_g2.dat');


%% TRACE FIG 3B
color_gray = [ 0.3 0.3 0.3]; 
T1=3870/0.01; 
T2= 4770/0.01;  

figure
hold on
plot(t(T1:T2), VI2(T1:T2,1),'color', color_gray, 'linewidth', 1)
    box off
    xticks(0)
    xticklabels({''})
    yticks([-50 0])
    yticklabels({'',''})
    ylim([-85 55])
    axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx pty]);
    print('/Users/kathleen/Documents/PhD/2022-reset/Fig3_trace/Fig3_VE_var','-depsc')%, '-painters')
%%
figure
hold on
plot(t(T1:T2), VI2(T1:T2,3),'color', color_gray, 'linewidth', 1)
    box off
    xticks(0)
    xticklabels({''})
    yticks([-50 0])
    yticklabels({'',''})
    ylim([-85 55])
    axis off
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 ptx pty]);
    print('/Users/kathleen/Documents/PhD/2022-reset/Fig3_trace/Fig3_VC_var','-depsc')%, '-painters')

    