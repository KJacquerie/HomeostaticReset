function function_fig_wf(f, data_wneg,data_wpos, name, ymax)
pt=11;
figure
hold on
plot([0 50], [1 1 ], '-', 'color', [0.9 0.9 0.9], 'linewidth', 0.5)
plot(f,data_wneg,':', 'color', 'k', 'linewidth', 1) 
plot(f,data_wpos,'-','color', 'k', 'linewidth', 1) 
%errorbar(Sjo(1:2:end,1),Sjo(1:2:end,3)+1, Sjo(1:2:end,4),'s','color', [0.5 0.5 0.5],'markersize',2.8)
%errorbar(Sjo(2:2:end,1),Sjo(2:2:end,3)+1, Sjo(2:2:end,4),'o','color', [0.5 0.5 0.5],'markersize',2.8)
%le= legend('75','','SJO','');

xlim([-1 51])
ylim([0.4 ymax])
%ylim([0.4 2])
%set(le, 'fontsize',pt, 'interpreter','latex', 'location', 'northwest')
%xlabel('$f [Hz]$', 'fontsize', pt, 'interpreter','latex')
%ylabel('$\Delta w [-]$', 'fontsize', pt, 'interpreter','latex')
set(gcf,'PaperPositionMode','auto');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 5 3.5]);
print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig_datasheet_wf/wf_%s',name), '-depsc');


end