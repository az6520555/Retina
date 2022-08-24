%%% test sustained/transient response with LN model, Rona, May 2018 %%%

%% stimulus
%ONOFF
st = -1*ones(1,6000);
st(2001:4000) = 1;
%HMM
Tot = 10000;  
dt = 50*10^-3; 
dtau = 50*10^-3; %step width 
T = dtau:dtau:Tot;
G = 10; % damping  
w = G/(2*1.06);%w = G/(2w)=1.06;  
D = 4; %dynamical range
a2=[];ey=[];ey1=[];
L = zeros(1,length(T));
V = zeros(1,length(T));
for t = 1:length(T)
    L(t+1) = L(t) + V(t)*dt;  %(th*(mu - X(t)))*dt + sqrt(dt)*randn; th = .01;  mu = 3;
    V(t+1) = (1-G*dt)*V(t) - w^2*L(t)*dt + sqrt(dt*D)*randn;
end      
st = L-mean(L);
st = st/(max(st)-min(st))*2;

% figure(2);plot(st);
%% linear filter
% STA = -1*gampdf(0:450,14,10)+1.4*gampdf(0:450,14,14);  % transient:-1.5,1.4,14,10,14,14
STA = -1.5*gampdf(0:450,15,12)+0.8*gampdf(0:450,15,15); % sustained:-1.5,0.8,15,12,15,15
% STA = STA./max(abs(STA));
figure(1);hold on;plot(STA)

F = conv(st,STA);
F = F(1:length(st));
F = F/max(F);
F(F<0.3) = 0;
figure(2);hold on;plot(F);
%% Nonlinear Function
g = 1;
th = 5;
rr = 100;
xx = 0:0.01:20;
G = 1./(1+rr*exp(-g*(xx-th)));
figure(4);plot(xx,G)
LN = 1./(1+rr*exp(-g*(F-th)));

figure(5);hold on;plot(LN)

%% MI
[t,mi] = MIfun(F,st);
figure(6);hold on;plot(t,mi/max(mi));

