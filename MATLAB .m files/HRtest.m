clear all
close all
clc

Tsmall = load('Tsmall');%zeros(2,8);  
Ttrans = 500; 

%wHR_Tsmall = zeros(1,8);
%wHR_Tfull = zeros(1,8);
%wHR_simu = zeros(1,8);
%wSAT_Tsmall = zeros(1,8);

%% HOMEOSTATIC RESET

idx_save=5; % choix du courant à tester

V = load(sprintf('V%d.dat',idx_save));
% w = load(sprintf('w%d.dat', idx_save)); 
% tu ne le loades pas encore car tu
% fais les test pour anticiper la valeur de HR

%% save value of w 
t=0.01:0.01:80000; 
% figure
% plot(w)

% wHR_simu(idx_save) = w(end); 
% à vérifier qu'on est en dehors du burst 
% (idx_save=4: 0.8104; idx_save = 7: 0.6648)

% save('wHR_simu.mat'); 


%%


test=0;  % 0= full, 1=small window


if(test==1)
T1=Tsmall(1,1);
T2= Tsmall(1,2);
else
    T1 = 10000/0.01;
    T2 = 80000/0.01; 
end
V1 = V(T1:T2,1);
V3 = V(T1:T2,3);

figure
hold on
plot(V1)
plot(V3)

%%

if(test==1)
    Tsmall(idx_save,1) = T1; 
    Tsmall(idx_save,2) = T2; 
end


Vspk1 = zeros(length(V1),1); 

Vspk3 = Vspk1; 

Vspk1(1)= 0; 
Vspk3(1)=0; 

for idx=2:1:length(V1)
    if(V1(idx)>0 && V1(idx-1)<0)
        Vspk1(idx) =1; 
    end
    if(V3(idx)>0 && V3(idx-1)<0)
        Vspk3(idx) =1; 
    end
    
end


[c, lags] = xcorr(Vspk3,Vspk1);
 
figure
stem(lags*0.01,c./(T2*0.01-T1*0.01))


A_p = 0.0096;
A_m = 0.0053; 
tau_p = 16.8;
tau_m = 33.7; 

A = zeros (2, length(lags));
A(1,:) = lags(:); % index
A(2,:) = c(:); % valeurs de la correlation 

F1=0; 
F2=0; 

% calcul de C+
for s = floor(length(A)/2)+1:1:length(A)
    
    F1 = F1 + exp(-((A(1,s))*0.01)/tau_p)*A(2,s);

end 
% calcul de C-
for s = 1:1:floor(length(A)/2)+1
    
    F2 = F2 + exp(((A(1,s))*0.01)/tau_m)*A(2,s);

end 
   


mu = 1;
lambda = (A_p/A_m)*(F1/F2);
wHR=lambda / (1+lambda);
wmu = (1+(1/(lambda^(1/mu))))^-1;
wSAT = (A_p*F1 - A_m*F2)/((T2-T1)*0.01);