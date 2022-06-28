%clear all
%close all
clc

pt=11;

%set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
%set(groot, 'defaultLegendInterpreter','latex');



color_gray = [0.3 0.3 0.3]; 

%%
dt=0.01;

V_wake = load('D:/PhD/Caro_Fig4/V_wake.mat');
V_sleep = load('D:/PhD/Caro_Fig4/V_sleep.mat');


%V_wake = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/%s/Results_PairBased/Vconnect_spk.dat','NetBig_wake'));

%V_sleep = load(sprintf('/Users/kathleen/Documents/PhD/2022-reset/%s/Results_PairBased/Vconnect_spk.dat','NetBig'));

%%
nPre = 100; 
nPost = 100; 
T1 = 10000/dt;
T2 = 11000/dt;

Vspk1 = zeros(T2-T1,nPre); % PRE
Vspk3 = Vspk1; %zeros(size(V_wake,2),1); %POST

idx_V = floor(V_wake.V_wake*1/dt);


%% Construction du vecteur 0 quand spike pas et 1 spike 

for i=1:1:nPre
    for curs=1:1:size(idx_V,2) 
        if(idx_V(i,curs) >=T1)
            if(idx_V(i,curs) <=T2 )
                
                Vspk1(idx_V(i,curs)-T1,i) = 1;
            end
        end
    end
end
%%
curs=1; 
for i=nPre+1:1:(nPre+nPost)
    for curs=1:1:size(idx_V,2) 
        if(idx_V(i,curs) >= T1 && idx_V(i,curs)<=T2 )
            Vspk3(idx_V(i,curs)-T1,i-nPre) = 1;
        end
    end
end

%%

scatter_V(Vspk1,Vspk3, 'Wake_200', 'large');

