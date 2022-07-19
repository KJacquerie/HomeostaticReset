function function_plot_wtimePrediction(folder, expm, strength)
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    f='/Users/JuliettePonnet/Documents/MASTER2/TFE/Code/Graupner2016/'  ;
%%
    w0 = 0:0.1:1; 
    dt=0.01;
    IstepI = -4.7:0.1:-4;
    T = 80000;
    n_loops = length(w0)*length(IstepI);
    
     %%

    wplot=zeros(n_loops, T);
    wf= zeros(1,n_loops); 
    if folder == 'GB_DD'
        demo = zeros(n_loops,4);
    else
       demo = zeros(n_loops,5) ;
    end
 
    
    idx=1;
    for m=1%:1:length(w0)
            for p=1%:1:length(IstepI)
            %wplot(idx,:) = load(sprintf('%s%s/results/results_%s/%s/w%d/%d/w.dat',f,folder,expm,strength, floor(w0(m)*10), p));
            %wf(idx)      = load(sprintf('%s%s/results/results_%s/%s/w%d/%d/wf.dat',f,folder,expm,strength, floor(w0(m)*10), p));
            demo(idx,:) = load(sprintf('%s%s/results/results_%s/BurstMATD/%s/w%d/%d/demo.dat',f,folder,expm,strength, floor(w0(m)*10), p));
            idx=idx+1;      
            end
    end
    %% Predict w using demo 
    
     w_predict = predictw(demo,expm); 
    %% prediction w (soft bound)
    color_test = {}; 
    color_test{1}=[0,  32,96,];
    color_test{2} =[42,  70, 122,]; 
    color_test{3}=[88, 110,150]; 
    color_test{4}=[ 109, 136,183]; 
    color_test{5}=[143,  170, 220,]; 
    color_test{6}=[ 180,  199,231,]; 
    color_test{7}=[ 205,  224,243,]; 
    color_test{8}=[ 223,  231,245,]; 

    figure
    for i=1:1:n_loops
        xx = wplot(i,:);
        hold on 
        plot(xx', 'color',  color_test{mod(i,length(IstepI))+1}/255,'linewidth',1.9)
        plot(length(xx)-10000,wf(i),'x','MarkerSize',10, 'MarkerEdgeColor', color_test{mod(i,length(IstepI))+1}/255)
        plot(length(xx)-10000,w_predict(i),'o','MarkerSize',10, 'MarkerEdgeColor', color_test{mod(i,length(IstepI))+1}/255)
        xlabel('t')
        ylabel('w')
    end
    
    xlabel('t', 'Interpreter', 'latex', 'fontsize', 11)
    ylabel('w', 'Interpreter', 'latex', 'fontsize', 11)
    title('')
   % ylim([0.5 0.8])
    ylim([0 1])
    %xlim([62000 72000])
    
    
    set(gcf, 'defaultAxesTickLabelInterpreter','latex'); 
    set(gcf, 'defaultLegendInterpreter','latex');
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 9.5 8]);
%    print(sprintf('/Users/JuliettePonnet/Documents/MASTER2/TFE/Code/Graupner2016/%s/results/figures_%s/BurstMATD/%s/wtime_pred',folder,expm,strength),'-depsc')

%%  prediction slope (hard bound)
    color_test = {};   
    color_test{1}=[0,  32,96,];
    color_test{2} =[42,  70, 122,]; 
    color_test{3}=[88, 110,150]; 
    color_test{4}=[ 109, 136,183]; 
    color_test{5}=[143,  170, 220,]; 
    color_test{6}=[ 180,  199,231,]; 
    color_test{7}=[ 205,  224,243,]; 
    color_test{8}=[ 223,  231,245,]; 
    
    color_pred = {};
    color_pred{1} = [236,236,236];
    color_pred{2} = [225,225,225];
    color_pred{3} = [215,215,215];
    color_pred{4} = [205,205,205];
    color_pred{5} = [195,195,195];
    color_pred{6} =[175,175,175];
    color_pred{7} =[145,145,145];
    color_pred{8} =[115,115,115];
    color_pred{8} = [100,100,100];

    figure
    
    for i=1:1:n_loops
        xx = wplot(i,:);
        hold on 
        plot(xx', 'color',  color_test{mod(i,length(IstepI))+1}/255,'linewidth',0.5)
        xlabel('t')
        ylabel('w')
    end
    
    hold on
 
    num_w = 0;
    idx_color=1;
    for i=length(IstepI)*num_w+1:1:(length(IstepI)*(num_w+1))
         hold on
         xx = wplot(i,:);
         xxx = 0:0.25:length(xx); %just a vector to plot the trace curve predicted
        if demo(i,5) %if at the end of the simulation w reached 1 
            x1 = demo(i,5)*dt; %coordinate of first point where w=1 to have the starting point for the curve to trace
            y1 = 1;
        else %if at the end of the simulation w did not yet reached 1 
            x1 = length(xx);
            y1 = xx(end);
        end
        m = w_predict(i);
        y = m*(xxx - x1) + y1; % curve with the slope predicted
        plot(xxx,y, 'color', color_pred{mod(idx_color,length(IstepI))+1}/255,'LineWidth',10)
        plot(xx', 'color',  color_test{mod(idx_color,length(IstepI))+1}/255,'linewidth',1.9)
        idx_color=idx_color+1;
    end
     ylim([0 1])

    
    set(gcf, 'defaultAxesTickLabelInterpreter','latex'); 
    set(gcf, 'defaultLegendInterpreter','latex');
    set(gcf,'PaperPositionMode','auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 9.5 8]);
%    print(sprintf('/Users/JuliettePonnet/Documents/MASTER2/TFE/Code/Graupner2016/%s/results/figures_%s/BurstMATD/%s/wtime_pred',folder,expm,strength),'-depsc')
end
