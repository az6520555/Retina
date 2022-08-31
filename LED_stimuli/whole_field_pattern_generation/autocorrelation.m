%finding autocorrelation of HMM and OU
clear
close all

L=[]
K=[]
% ============ HMM ==============
Tot = 300;  
dt = 50*10^-3; 
dtau = 50*10^-3; %step width 
T = dtau:dtau:Tot;
G = 5; % damping  
w = G/(2*1.06);%w = G/(2w)=1.06;  
D = 4; %dynamical range
at = 3;
m = 0.05;

a2=[];ey=[];ey1=[];
L = zeros(1,length(T));
V = zeros(1,length(T));
for t = 1:length(T)
    L(t+1) = L(t) + V(t)*dt;  %(th*(mu - X(t)))*dt + sqrt(dt)*randn; th = .01;  mu = 3;
    V(t+1) = (1-G*dt)*V(t) - w^2*L(t)*dt + sqrt(dt*D)*randn;
end    
plot(T,L(1:end-1))
% ===============================

% ============= OU ===============
Tot = 300;  
dt = 50*10^-3; 
dtau=50*10^-3; %step width 
T = dtau:dtau:Tot;
tau=1;
D = 4; %dynamical range
at = 30;
m = 0.05;
a2=[];ey=[];ey1=[];
K = zeros(1,length(T));
for i = 1:length(T)-1
K(i+1) = K(i) + (-K(i)/tau + randn*sqrt(D/dt))*dt;
end
figure;plot(T,K)


[ACF1,lags1]=autocorr(L);
[ACF2,lags2]=autocorr(K);
% HMM correlation time
[M ind1]=min(abs(ACF1-0.5));
corrt1=ind1/dtau
% OU correlation time
[M ind2]=min(abs(ACF2-0.5));
corrt2=ind2/dtau
