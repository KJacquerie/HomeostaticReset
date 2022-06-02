function function_fig_STDP(dt, f, data, name)

    switch f
        case 1
            color_vec = [0.5 0.5 0.5];
        case 50
            color_vec = [0.1 0.1 0.1];
    end
 

    pt=11;
    figure
    hold on
    plot([-100 100], [1 1], '-', 'color', [0.8 0.8 0.8], 'linewidth',0.5)
    plot([0 0], [0 2], '-', 'color', [0.8 0.8 0.8], 'linewidth',0.5)
    plot(dt, data, 'color', color_vec)
    
    xlim([-80 80])
    %ylim([0.4 1.6])
    %xlabel('$\Delta t [ms]$', 'fontsize', pt, 'interpreter','latex')
    %ylabel('$\Delta w [-]$', 'fontsize', pt, 'interpreter','latex')
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5 3.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig_datasheet_STDP/STDP_%s_f%d',name,f), '-depsc');


end
