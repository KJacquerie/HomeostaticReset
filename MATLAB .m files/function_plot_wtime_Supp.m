function function_plot_wtime_Supp(folder, expm, name_save,strength)
    w0 = 0:0.1:1; % 11 items
    
    IstepI = -4.7:0.1:-4; % 8 items

    w0_mat = zeros(length(w0), length(IstepI));
    val = 0;
     for i=1:1:length(w0)
          w0_mat(i,:) = val ; 
          val=val+0.1;
     end
     %%
    wplot=zeros(88, 80000); 
    t=1:1:80000;
    wf= zeros(1,88); 
    wf_mat = zeros(length(w0), length(IstepI)); 
    
    idx=1;
    for m=1:1:length(w0)
        for p=1:1:length(IstepI)
            wplot(idx,:) = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/%s/results_%s/%s/w%d/%d/w.dat',folder,'sleep',strength, floor(w0(m)*10), p));
            wf(idx)      = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/%s/results_%s/%s/w%d/%d/wf.dat',folder,'sleep',strength, floor(w0(m)*10), p));
            wf_mat(m,p)  = wf(idx); 
            idx=idx+1;

        end
    end

    %%

    color_test = {}; 
    color_test{1}=[8,  37,103,];
    color_test{2} =[0,  0, 128,]; 
    color_test{3}=[0, 24,168]; 
    color_test{4}=[ 0,  0,255,]; 
    color_test{5}=[ 0, 112,187]; 
    color_test{6}=[118,  171, 223,]; 
    color_test{7}=[ 70,  102,255,]; 
    color_test{8}=[ 108,  191,238,]; 
    %color_test{9}=[  175,  219,245]; 


    figure
    for i=1:1:88
        xx = wplot(i,:);
        hold on 
        plot(t*1e-3,xx', 'color',  color_test{mod(i,8)+1}/255)
        %xlabel('t')
        %ylabel('w')
    end

    %xlabel('t', 'Interpreter', 'latex', 'fontsize', 14)
   % ylabel('w', 'Interpreter', 'latex', 'fontsize', 14)
    title('')

    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 5 3.5]);
    print(sprintf('/Users/kathleen/Documents/PhD/2022-reset/Fig_datasheet_reset/wtime_%s_%s_supp',name_save, strength),'-depsc', '-painters')
    %print(sprintf('wtime_%s',expm),'-dpdf')
    
    
    %%
%     color_test = {}; 
%     color_test{1} =[255,  230, 153,]; 
%     color_test{2} =[238,  185, 115,]; 
%     color_test{3}=[232, 166,79]; 
%     color_test{4}=[ 211,  134,50,]; 
%     color_test{5}=[ 195,  212,214,]; 
%     color_test{6}=[133,  168,172,];
%     color_test{7}=[ 86, 127,129]; 
%     color_test{8}=[ 69,  113,115]; 
%     color_test{9}=[ 46,  82, 83]; 
%     color_test{10}=[  33,  59,60]; 
%        color_test{11}=[  0,  0,0]; 
%     %%
%     figure
%     idx_color=1;
%     for i=6:9:91
%         xx = wplot(i,1:40000);
%         hold on 
%         plot(xx', 'color',  color_test{idx_color}./255)
%         xlabel('t')
%         ylabel('w')
%         idx_color=idx_color+1;
%     end
% 
%     xlabel('t', 'Interpreter', 'latex', 'fontsize', 14)
%     ylabel('w', 'Interpreter', 'latex', 'fontsize', 14)
%     title('')
% 
%     set(gcf,'PaperPositionMode','auto');
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf, 'PaperPosition', [0 0 18 15]);
%     print(sprintf('wtime_%s','DOneline'),'-depsc')
%     print(sprintf('wtime_%s','DOneline'),'-dpdf')
end
