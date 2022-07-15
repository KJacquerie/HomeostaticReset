
function scatter_simuvspredict(expm_state, Twindow)
    color_blue = [106 153 208]./255; 
    nPre = 100; 
    nPost = 100; 
    y = load(sprintf('D:/PhD/Caro_Fig4/%s_large_wsimu.mat', expm_state));
    
    wsimu = reshape(y.wsimu, nPre*nPost,1); 

    y = load(sprintf('D:/PhD/Caro_Fig4/%s_%s_wpredict.mat', expm_state, Twindow));
    wpredict = reshape(y.wpredict, nPre*nPost,1); 



    %%


    figure

    hold on
    scatter(wsimu, wpredict, 1, 'MarkerEdgeColor',color_blue,'MarkerFaceColor',color_blue)
    xlim([min(wsimu) max(wsimu)])
    ylim([min(wpredict) max(wpredict)])
    xticks([0 1])
    xticklabels({'',''})
    yticks([0 1])
    yticklabels({'',''})
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 1.5 1]);
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/scatterZoom_%s_%s',expm_state, Twindow ), '-dmeta', '-painters')
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/scatterZoom_%s_%s',expm_state, Twindow ), '-depsc', '-painters')


    hold on
    plot([0 1], [0 1],'color', [0.5 0.5 0.5])

    xlim([0 1])
    ylim([0 1])
    %xlabel('w simu')
    %ylabel('w predict')
    %axis off
    
    xticks([0 1])
    xticklabels({'',''})
    yticks([0 1])
    yticklabels({'',''})
    
    xlim([0 1])
    ylim([0 1])
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 2]);
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/scatter_%s_%s',expm_state, Twindow ), '-dmeta', '-painters')
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/scatter_%s_%s',expm_state, Twindow ), '-depsc', '-painters')


    
    

end
