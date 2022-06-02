clear all
close all 
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

ptx = 2.25;
pty = 0.75;

for idx_trace=5:1:5
   
    
t=0.01:0.01:50000; 

w = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Network_heterogeneous/trace%d/w.dat',idx_trace));
V = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Network_heterogeneous/trace%d/Vconnect.dat',idx_trace));
%%

color_w={}; 
color_w{1} = [44, 134, 237]; 
color_w{2} = [99, 205, 191] ; 
color_w{3} = [219, 109, 139]; 
color_w{4} = [241, 201, 103]; 
color_w{5} = [205, 122, 62]; 
color_w{6} = [156, 63, 140]; 
color_w{7} = [66, 148, 135]; 
color_w{8} = [39, 82, 177]; 
color_w{9} = [77, 10, 134]; 

%%
T1=31100/0.1;
T2= 39100/0.1;
color_gray = [ 0.3 0.3 0.3]; 
for idx_V=1:1:7
    figure
    plot(V(T1:T2,idx_V), 'color', color_gray, 'linewidth', 1)
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
    print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig3_ntk/ntk%d_V%d',idx_trace, idx_V),'-depsc')%, '-painters')

end
%%

figure
hold on 
for idx_color=1:1:size(w,2)
    plot(t,w(1:end,idx_color),'color',  color_w{idx_color}./255, 'linewidth',1)
end
xticks([0 1 2 3 4 5]*10000) 
xticklabels({'','','','','',''})
yticks([0 0.5 1]) 
yticklabels({'','',''})
%xlabel('t', 'Interpreter', 'latex', 'fontsize', 14)
%ylabel('w', 'Interpreter', 'latex', 'fontsize', 14)
title('')
ylim([-0.05 1.05])
xlim([-10 t(end)+10])

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6 3.2]);

print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig3_ntk/w%d',idx_trace),'-depsc', '-painters')
%print(sprintf('wtime_%s',expm),'-dpdf')

end