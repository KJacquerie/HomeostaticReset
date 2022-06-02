clear all   
close all
clc
folder = "model2" ;
expm = "Control" ;
strength = "weak";
%%
w0 = 0:0.1:1; 

IstepI = -4.7:0.1:-4;

n_net = 10;
T = 80000;
n_loops = length(w0)*length(IstepI);
%%

wplot=zeros(n_loops, T,n_net); 

wf= zeros(n_loops,n_net); 

idx=1;
for m=1:1:length(w0)
    for p=1:1:length(IstepI)

        wplot(idx,:,:) = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model2/netBurstMATD/%s/w%d/%d/w_net.dat',strength, floor(w0(m)*10), p));
        wf(idx,:)      = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/model2/netBurstMATD/%s/w%d/%d/wf_net.dat',strength, floor(w0(m)*10), p));

        idx=idx+1;
    end
end

wf_reshaped=reshape(wf,[length(IstepI),length(w0),10]);

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


color_pastel = {}; 
color_pastel{1} =[244,  229, 163,]; 
color_pastel{2} =[247,  208, 133,]; 
color_pastel{3}=[237, 170,99];  
color_pastel{4}=[ 229,  135,117,]; 
color_pastel{5}=[228,  135,148,];
color_pastel{6}=[ 230, 150,185]; 
color_pastel{7}=[ 221,  185,211]; 
color_pastel{8}=[ 138,  139, 191]; 
color_pastel{9}=[ 159 ,  201,224]; 
color_pastel{10}=[  172,  215, 216]; 
color_pastel{11}=[  163,  207, 172]; 
color_pastel{12}=[  205,  224, 147];

%%

idx_i = 8;
wI = wplot(idx_i:8:end, :,:); 

mean_w = mean(wI(:,:,:),3); 


std_w = std(wI(:,:,:),1,3);

%
figure
for idx_w =1:1:11
hold on
    errorbar(mean_w(idx_w,1:40000)', std_w(idx_w,1:40000)', 'color', color_pastel{idx_w}./255); 
end
plot(mean_w(:,1:40000)','k', 'linewidth', 1.5); 
ylim([0 1])

%%

expm = 'rainbow'; 
switch expm 
    case 'rainbow'
        name_save='_VARIABILITY_rainbow'; 
        color_shadow = color_pastel; 
    case 'blue'
        name_save='_VARIABILITY_blue'; 
        color_shadow = color_blue; 
end
figure
t=1:1:40000;
for idx_w=1:1:11
    curve1 = (mean_w(idx_w, 1:40000)'+ std_w(idx_w,1:40000)')';
    curve2 = (mean_w(idx_w, 1:40000)'- std_w(idx_w,1:40000)')';
    
hold on
box off
fill([t fliplr(t)], [curve1  fliplr(curve2)], color_shadow{idx_w}./255, 'EdgeColor', 'none', 'FaceAlpha', 0.5);

end
for idx_w=1:1:11
    plot(mean_w(idx_w,1:40000)','color','k' , 'linewidth', 1); 
end
ylim([0 1])

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
    
    
    
