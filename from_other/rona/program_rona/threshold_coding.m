% Simulate N and P cell by threshold detection, by Rona, 2018 
clearvars -except noise
rate = 1000;%HZ
Tot = 1000;%s  
dt = 0.05;%s 
dtau = 0.05; %step width 
T = dtau:dtau:Tot;
m = 0.05;
D = 4; %dynamical range

Stim = input('HMM(h)/OU(o)?  ')

if Stim == 'h'
    G = 0.5; % damping  
    w = G/(2*1.06);%w = G/(2w)=1.06;  
    a2=[];ey=[];ey1=[];
    L = zeros(1,length(T));
    V = zeros(1,length(T));
    for i = 2:length(T)
        L(i+1) = L(i) + V(i)*dt; 
        V(i+1) = (1-G*dt)*V(i) - w^2*L(i)*dt + sqrt(dt*D)*noise(i);
    end  

elseif Stim == 'o'
    tau = 7.5;
    a2=[];ey=[];ey1=[];
    L = zeros(1,length(T));
    for i = 2:length(T)
        L(i+1) = L(i) + (-L(i)/tau + noise(i)*sqrt(D/dt))*dt;
    end 
    L = smooth(L,50)';
end

L = 0.01*L/std(L);
X=L;% X = L-mean(L)+m; % X is the final stimuX
figure;autocorr(X,500)
figure(10);hold on;plot(X)

%% threshold detection
m = mean(X);
d = std(X);
spike = zeros(1,length(X));
for i = 1:length(X)-1
    if  m+0.8*d<X(i) && X(i)< m+1.2*d %N(on) cell
%     if X(i+1)-X(i)<0 && m-1.2*d<X(i) && X(i)<m-0.8*d %P(off) cell
        spike(i) = 1;
    end
end

%% MI
S = X;
figure;autocorr(S,100)
figure;plot(S/max(S));hold on;plot(spike)
[t,information,corr] = MIfun(spike,S,1,200,200);
figure(100);hold on;plot(t,information);%plot(locs,pks,'*')

