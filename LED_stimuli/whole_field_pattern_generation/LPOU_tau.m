% test OU off different tau
clear
close all
Tot = 300;  
dt = 0.01; 
T = dt:dt:Tot;
tau=[0.1 0.5];
D =4; %dynamical range
at = 60;
m = 10;
a2=[];ey=[];ey1=[];

random_term=randn(1,length(T))*sqrt(D/dt)*dt;
L1 = zeros(1,length(T));
for i = 1:length(T)-1
    L1(i+1) = L1(i) + (-L1(i)/tau(1) + random_term(i))*dt;
end  
L2 = zeros(1,length(T));
for i = 1:length(T)-1
    L2(i+1) = L2(i) + (-L2(i)/tau(2) + random_term(i))*dt;
end
figure(1);hold on;plot(T,L1);
figure(1);plot(T,L2);
legend(['tau=',num2str(tau(1))],['tau=',num2str(tau(2))])
xlabel('time (s)')

    %% fft
Fs = 1/dt
Ts =  dt            % Sampling period       
Len = length(L1);             % Length of signal
t = (0:Len-1)*Ts;        % Time vector
Y = fft(L1);
P2 = abs(Y/Len);
P1 = P2(1:Len/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(Len/2))/Len;
amp1=P1;
figure(2);hold on;plot(f,amp1) 
xlabel('f (hz)')
ylabel('power')
xlim([0 5])

Fs = 1/dt
Ts =  dt            % Sampling period       
Len = length(L2);             % Length of signal
t = (0:Len-1)*Ts;        % Time vector
Y = fft(L2);
P2 = abs(Y/Len);
P1 = P2(1:Len/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(Len/2))/Len;
amp2=P1;
figure(2);hold on;plot(f,amp2) 
legend(['tau=',num2str(tau(1))],['tau=',num2str(tau(2))])
xlim([0 5])

gain=amp2./amp1;
figure(3);plot(f,gain)
xlim([0 5])
xlabel('f')
ylabel('times')
