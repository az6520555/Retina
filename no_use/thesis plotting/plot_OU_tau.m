clear
close all
Tot = 300;  
dt = 0.01; 
T = dt:dt:Tot;
tau1=0.1;
tau2=1;
D =4; %dynamical range
at = 60;
a2=[];ey=[];ey1=[];
saving=0;

randseed=randn(1,length(T));
L1 = zeros(1,length(T));
L2 = zeros(1,length(T));
for i = 1:length(T)-1
    L1(i+1) = L1(i) + (-L1(i)/tau1 + randseed(i)*sqrt(D/dt))*dt;
    L2(i+1) = L2(i) + (-L2(i)/tau2 + randseed(i)*sqrt(D/dt))*dt;
end
[val1,lags1] =autocorr(L1,300);
[val2,lags2] =autocorr(L2,300);

plot(T,L1/std(L1),'linewidth',1);hold on;
plot(T,L2/std(L2),'linewidth',1)
legend('\tau=0.1s','\tau=1s')

figure;hold on;box on
plot(lags1,val1,'linewidth',1)
plot(lags2,val2,'linewidth',1)
legend('\tau=0.1s','\tau=1s')
xlabel('lags')
ylabel('autocorrelation')

figure
subplot(1,7,1:5);box on;hold on;
plot(T,L1/std(L1),'linewidth',1)
plot(T,L2/std(L2),'linewidth',1)
legend('\tau=0.1s','\tau=1s')
xlabel('time (s)')
subplot(1,7,6:7);box on;hold on;
plot(lags1,val1,'linewidth',1)
plot(lags2,val2,'linewidth',1)
xlabel('lags')
ylabel('autocorrelation')

Fs = 1/dt;
Ts =  dt       ;     % Sampling period       
Len = length(L2);             % Length of signal
t = (0:Len-1)*Ts;        % Time vector
Y2 = fft(L2);
P2 = abs(Y2/Len);
P1 = P2(1:Len/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(Len/2))/Len;
figure(5);hold on;box on
plot(f,P1,'linewidth',1,'color',[0.4660, 0.6740, 0.1880]);

Len = length(L1);             % Length of signal
t = (0:Len-1)*Ts;        % Time vector
Y1 = fft(L1);
P2 = abs(Y1/Len);
P1 = P2(1:Len/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(Len/2))/Len;  
plot(f,P1,'linewidth',1,'color',[0.6350, 0.0780, 0.1840]);
title('FFT')
xlabel('f (Hz)');ylabel('Amplitude')
xlim([0 5])
set(gcf,'Position',[300,300,300,300])
legend('\tau=1s','\tau=0.1s')