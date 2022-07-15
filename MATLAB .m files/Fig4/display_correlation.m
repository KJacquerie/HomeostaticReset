function display_correlation(expm_state, Twindow)
    
    y = load(sprintf('D:/PhD/Caro_Fig4/%s_%s_corr_neg.mat', expm_state, Twindow));
    Corr_neg = y.Corr_neg; 
    y = load(sprintf('D:/PhD/Caro_Fig4/%s_%s_corr_pos.mat', expm_state, Twindow));
    Corr_pos = y.Corr_pos; 

    %% Color 
    
 

%     newmap = jet;                    %starting map
%     ncol = size(newmap,1);           %how big is it?
%     zpos = 1 + floor(2/3 * ncol);    %2/3 of way through
%     newmap(zpos,:) = [1 1 1];        %set that position to white
%     colormap(newmap);                %activate it

    %%
    Corr = Corr_pos./Corr_neg;

  
    figure; imagesc(Corr); colormap( [1 1 1; parula(256)] );caxis([0 1.7]);
    %colorbar;
    %set(gca,'ColorScale','log')
    %caxis([log10(0.00001) log10(1e13)]);%caxis([0 1])

    %xlabel('$\#$ pre', 'interpreter', 'latex')
    %ylabel('$\#$ post', 'interpreter', 'latex')
    
    xticks([0 50 100])
    xticklabels({'','', ''})
    yticks([0 50 100])
    yticklabels({'','',''})
    
    xlim([0 100])
    ylim([0 100])
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 2]);
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/corr_%s',expm_state ), '-dmeta', '-painters')
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/corr_%s',expm_state ), '-depsc', '-painters')
    
    figure
    colorbar;
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 2]);
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/bar_corr_%s',expm_state ), '-dmeta', '-painters')

end