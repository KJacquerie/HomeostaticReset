%% HOMEOSTATIC RESET

%T1=zeros(1,8);
%T2=zeros(1,8);

%%
idx_save=8; 

w = load(sprintf('wHB%d.dat', idx_save));

plot(w)
wmat(idx_save,:)=w; 

%T1(idx_save) = 732;
%T2(idx_save) = 11810;

wSAT_test(idx_save) =( w(T2(idx_save))-w(T1(idx_save)))/(T2(idx_save)- T1(idx_save)); 


%%
%T1bis =3000* ones(1,8); 
%wSAT_bis = ( w(T2(idx_save))-w(T1bis(idx_save)))./(T2-T1bis);
%%
% for i = 1:1:8
% wm(i,:) = load(sprintf('wHB%d.dat', i))';
% 
% end

%%
tt= 3000:1:10000; 
wanal = [8.9977,5.0690,6.7638, 8.1860, 9.3400, 4.8326, 5.4758, -6.2585 ]; 

%%
figure
hold on


ii=3;
for ii=2:1:8
    vec = [T1(ii) T2(ii)];
plot(vec, wSAT_test(ii)*(vec-T1(ii))+wm(ii,T1(ii)), 'linewidth', 3)
end
plot(wm(1:8,:)')

figure
hold on

for ii=2:1:8
    vec = [T1(ii) T2(ii)];
plot(vec, wanal(ii)*1e-5*(vec-3000)+wm(ii,3000), 'linewidth', 3)
end
plot(wm(1:8,:)')
ylim([0 1]);
