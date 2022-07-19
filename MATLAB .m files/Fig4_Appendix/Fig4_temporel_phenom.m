close all
clear all;

pt=11; % taille de la police de votre rapport

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

strength = 'weak';
x = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/PairBased/trace/w_%s/xBurst.dat',strength));
y =  load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/PairBased/trace/w_%s/yBurst.dat',strength));
w =  load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/PairBased/trace/w_%s/wBurst.dat',strength));

%%
t=0.01:0.01:50000;

fig = 'T1';

if fig == 'T1'
    t1 = 1763/0.01;
    t2 = 2363/0.01;
else
    t1 = 49400/0.01;
    t2 = 50000/0.01;
end

figure;
plot(t(t1:t2),x(t1:t2),'Color',[0.8 0.8 0.8]);
hold on
plot(t(t1:t2),y(t1:t2),'Color',[0.5 0.5 0.5]);
axis off
ylim([0 12e-3])

set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5.5 1.5]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig4/PairBased_trace_%s_%s',strength,fig),'-depsc')

switch strength
    case 'weak' 
        switch fig
            case 'T1'
                y_vec=[0.115 0.14];
            case 'T2'
                y_vec=[0.66 0.685];
        end
    case 'strong'
        switch fig
            case 'T1'
                y_vec=[0.855 0.88];
            case 'T2'
                y_vec=[0.665 0.69];
        end
end

figure;
switch strength
    case 'weak'
    plot(t(t1:t2),w(t1:t2),'Color',[168, 212, 255]./255);
    axis off
    ylim(y_vec)
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5.5 1.5]);

    print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig4/PairBased_w_%s_%s',strength, fig),'-depsc')
    case 'strong'
    plot(t(t1:t2),w(t1:t2),'Color',[77, 135, 192]./255);
    axis off
    ylim(y_vec)
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5.5 1.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig4/PairBased_w_%s_%s',strength, fig),'-depsc')
end


% 
% 
% % Create subplot
% subplot1 = subplot(7,1,1 );
% plot(t,V1, 'Color',[0.474509805440903 0.725490212440491 0.976470589637756]);
% ylabel('E [$mV$]','FontSize',pt,'Interpreter','latex');
% box off
% hold on 
% subplot2 = subplot(7,1,2);
% plot(t,V2, 'Color', [1 0.400000005960464 0.400000005960464]);
% ylabel('$C$ [$mV$]','FontSize',pt,'Interpreter','latex');
% box off
% set(subplot1,'XTick',zeros(1,0));
% set(subplot2,'XTick',zeros(1,0));
% hold on 
% 
% subplot3 = subplot(7,1,3 );
% plot(t,x, 'Color',[0.474509805440903 0.725490212440491 0.976470589637756]);
% box off
% set(subplot3,'XTick',zeros(1,0));
% ylabel('$r_1$ [-]','FontSize',pt,'Interpreter','latex');
% hold on 
% 
% 
% subplot4 = subplot(7,1,4);
% plot(t,x2, 'Color',[0.474509805440903 0.725490212440491 0.976470589637756]);
% box off
% set(subplot4,'XTick',zeros(1,0));
% ylabel('$r_2$ [-]','FontSize',pt,'Interpreter','latex');
% hold on 
% 
% 
% 
% subplot5 = subplot(7,1,5);
% plot(t,y, 'Color',[1 0.400000005960464 0.400000005960464]);
% box off
% set(subplot5,'XTick',zeros(1,0));
% ylabel('$\sigma_1$ [-]','FontSize',pt,'Interpreter','latex');
% hold on 
% 
% subplot6 = subplot(7,1,6);
% plot(t,y2, 'Color',[1 0.400000005960464 0.400000005960464]);
% box off
% set(subplot6,'XTick',zeros(1,0));
% ylabel('$\sigma_2$ [-]','FontSize',pt,'Interpreter','latex');
% hold on 
% 
% subplot7 = subplot(7,1,7);
% plot(t,w, 'Color',[228, 184, 101]./255);
% box off
% ylabel('$w_{ij}$ [-]','FontSize',pt,'Interpreter','latex');
% xlabel('Time [ms]','FontSize',pt,'Interpreter','latex');
% hold on 
% 
% 
% set(gcf,'PaperPositionMode','auto');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 ptx pty]);
% % enregistre une figure avec un nom parametre
% print('Triplet_temp','-dpng')
% print('Triplet_temp','-depsc')
% print('Triplet_temp','-dmeta')
% 
% 
