clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');



%%
for idx_Ivec=1:1:9
        plot_BurstI(idx_Ivec)
end

function plot_BurstI (idx_Ivec)
V = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model2/trace/burst/Vconnect_%d.dat',idx_Ivec));

t1 = 0.01; 
t2 = 5000;
t=0.01:0.01:t2;
T1=2000/0.01; 
T2= 5000/0.01;  

figure
subplot(3,1,1)
plot(t(T1:T2), V(T1:T2,1),'color', [0.3 0.3 0.3], 'linewidth', 0.5)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
box off
ylim([-85 55])
%ylabel('$V_E$ [mV]', 'interpreter', 'latex', 'fontsize', 11)
yticks([-50 0 50])
yticklabels({'', '',''})
xticks([2000 3000 4000 5000])
xticklabels({'','','',''})


subplot(3,1,2)
plot(t(T1:T2), V(T1:T2,3),'color', [0.3 0.3 0.3], 'linewidth', 0.5)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
box off
ylim([-85 55])
%ylabel('$V_C$ [mV]', 'interpreter', 'latex', 'fontsize', 11)
yticks([-50 0 50])
yticklabels({'', '',''})
xticks([2000 3000 4000 5000])
xticklabels({'','','',''})


subplot(3,1,3)
plot(t(T1:T2), V(T1:T2,2),'color', [0.3 0.3 0.3], 'linewidth', 0.5)%, 'color', [168/255 12/255 64/255], 'linestyle', ':')
ylim([-85 55])
box off
yticks([-50 0 50])
yticklabels({'-50', '','50'})
xticks([2000 3000 4000 5000])
xticklabels({'0','1000','',''})
%ylabel('$V_I$ [mV]', 'interpreter', 'latex', 'fontsize', 11)
%xlabel('t [ms]', 'interpreter', 'latex', 'fontsize', 11)
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 4]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig_SI/Burst%d',idx_Ivec),'-depsc')%, '-painters')
end
