function function_plot_wtime_Article(folder, expm, name_save, strength, num_Iapp)
    w0 = 0:0.1:1; 

    IstepI = -4.7:0.1:-4;

    w0_mat = zeros(length(w0), length(IstepI));
    val = 0;
     for i=1:1:length(w0)
          w0_mat(i,:) = val ; 
          val=val+0.1;
     end
     %%

    wplot=zeros(88, 80000); 

    wf= zeros(1,88); 
    wf_mat = zeros(length(w0), length(IstepI)); 

    idx=1;
    for m=1:1:length(w0)
        for p=1:1:length(IstepI)
            wplot(idx,:) = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/%s/results_%s/%s/w%d/%d/w.dat',folder,expm,strength, floor(w0(m)*10), p));
            wf(idx)      = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/%s/results_%s/%s/w%d/%d/wf.dat',folder,expm,strength, floor(w0(m)*10), p));
            wf_mat(m,p)  = wf(idx); 
            idx=idx+1;

        end
    end

%%

color_blue = {}; 
color_blue{1} =[239,  243, 247,]; 
color_blue{2} =[224,  233, 240,]; 
color_blue{3}=[210, 222,233]; 
%color_blue{4}=[ 195,  211,226,]; 
color_blue{4}=[ 180,  199,219,]; 
color_blue{5}=[166,  188,211,];
color_blue{6}=[ 148, 174,203]; 
color_blue{7}=[ 137,  166,197]; 
color_blue{8}=[ 121,  153, 189]; 
color_blue{9}=[ 98 ,  134,177]; 
color_blue{10}=[  63,  106, 159]; 
color_blue{11}=[  56,  82, 121]; 


     color_test_Iapp = {}; 
    color_test_Iapp{1} =[240,  240, 240,]; 
    color_test_Iapp{2}=[220, 220,220]; 
    color_test_Iapp{3}=[ 190,  190,190,]; 
    color_test_Iapp{4}=[ 160,  160,160,]; 
    color_test_Iapp{5}=[120,  120,120,];
    color_test_Iapp{6}=[ 80, 80,80]; 
    color_test_Iapp{7}=[ 50,  50, 50]; 
    color_test_Iapp{8}=[ 10,  10, 10]; 

    
    %%
    figure(1)
    
    for i=1:1:88
        xx = wplot(i,1:40000);
        hold on 
        plot(xx', 'color',color_test_Iapp{mod(i,8)+1}/255, 'linewidth',0.25) 

        %if(mod(i,9)==0)
        %    idx_color = idx_color+1;
        %end
        
    end
    idx_color=1;
    for i=num_Iapp:8:88
        xx = wplot(i,1:40000);
        hold on 
        plot(xx', 'color',  color_blue{idx_color}./255, 'linewidth',1.5)
        %xlabel('t')
        %ylabel('w')
        idx_color=idx_color+1;
    end
    
    xticks([0 1 2 3 4]*10000) 
    xticklabels({'','','','',''})
    yticks([0 0.5 1]) 
    yticklabels({'','',''})
    %xlabel('t', 'Interpreter', 'latex', 'fontsize', 14)
    %ylabel('w', 'Interpreter', 'latex', 'fontsize', 14)
    title('')
    ylim([-0.05 1.05])

       set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 4 3.2]);
    
    print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig3_wtime/%s',name_save),'-depsc', '-painters')
    %print(sprintf('wtime_%s',expm),'-dpdf')
end

