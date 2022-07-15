function scatter_V(Vspk, Vspk3, expm_state, Twindow)
    color_blue = [106 153 208]./255; 
    nPre = 100; 
    nPost = 100; 




	%% figure
    
 for i=1:1:nPre
     hold on
     idx_plot = find(Vspk3(:,i)~=0);  
     for j=1:1:length(idx_plot)
         plot(idx_plot(j),i,'o', 'MarkerEdgeColor','k','MarkerFaceColor','k', 'Markersize', 0.05)        
     end  
 end
 
 for i=nPre+1:1:nPost+nPre
     hold on
     idx_plot = find(Vspk(:,i-nPre)~=0);  
     for j=1:1:length(idx_plot)
         plot(idx_plot(j),i,'o', 'MarkerEdgeColor','k','MarkerFaceColor','k', 'Markersize', 0.05)        
     end  
 end


    ylim([0 200])

    xticks([0 1])
    xticklabels({'',''})
    yticks([0 1])
    yticklabels({'',''})
    

    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 3 2]);
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/scatterV_%s',expm_state ), '-dmeta', '-painters')
    print(sprintf('/Users/carol.000/Documents/PhD/Data/fig/fig4/scatterV_%s',expm_state ), '-depsc', '-painters')


end