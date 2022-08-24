clear
close all
Tot = 300;  
dt = 0.01; 
T = dt:dt:Tot;
% tau=0.5;
D =4; %dynamical range
at = 60;
a2=[];ey=[];ey1=[];
saving=0;
tau1=0.5;
tau2=1.8;
tau3=20;
c={'k','r'}
L = zeros(1,length(T));
randseq=randn(1,length(T));
for i = 1:length(T)-1
    L(i+1) = L(i) + (-L(i)/tau1 + randseq(i)*sqrt(D/dt))*dt;
end
L2 = zeros(1,length(T));
for i = 1:length(T)-1
    L2(i+1) = L2(i) + (-L2(i)/tau2 + randseq(i)*sqrt(D/dt))*dt;
end
L2 = zeros(1,length(T));
for i = 1:length(T)-1
    L2(i+1) = L2(i) + (-L2(i)/tau2 + randseq(i)*sqrt(D/dt))*dt;
end
L3 = zeros(1,length(T));
for i = 1:length(T)-1
    L3(i+1) = L3(i) + (-L3(i)/tau3 + randseq(i)*sqrt(D/dt))*dt;
end


cutoff=1;
corrtimef=[];

[b,a]=butter(2,2*cutoff*dt,'low')
Lf=filter(b,a,L);

[val,lags] =xcorr(L,Lf,300);
[M,I]=max(val);
Lfshift=Lf(-lags(I)+1:end);

%normalization
L=L/std(L);
Lf=Lf/std(Lf);
L2=L2/std(L2);
L3=L3/std(L3);
Lfshift=Lfshift/std(Lfshift);

figure(1)
subplot(2,1,1)
plot(T,L,'linewidth',1);hold on;plot(T(1:end+lags(I)),Lfshift,'linewidth',1)
xlim([100,130])
% title('Normalized OU and LPOU')
legend('OU','LPOU, f_c=1Hz')
ylabel('$x(t)$','Interpreter','Latex')
% ylim([-3,3])
subplot(2,1,2)
plot(T(1:end-1),diff(L)/dt,'linewidth',1);hold on ;plot(T(1:end-1),diff(Lf)/dt,'linewidth',1)
xlim([100,130])
ylim([-80,80])
ylabel('$dx(t)/dt$','Interpreter','Latex')
xlabel('t (s)')

figure(2)
subplot(2,1,1)
plot(T,L,'linewidth',1);hold on;plot(T,L2,'linewidth',1);plot(T,L3,'linewidth',1)
xlim([100,130])
% title('Normalized OU and LPOU')
legend(['OU, \tau=',num2str(tau1),'s'],['OU, \tau=',num2str(tau2),'s'],['OU, \tau=',num2str(tau3),'s'])
ylabel('$x(t)$','Interpreter','Latex')
% ylim([-3,3])
subplot(2,1,2)
plot(T(1:end-1),diff(L)/dt,'linewidth',1);hold on ;plot(T(1:end-1),diff(L2)/dt,'linewidth',1);plot(T(1:end-1),diff(L3)/dt,'linewidth',1)
xlim([100,130])
ylim([-80,80])
ylabel('$dx(t)/dt$','Interpreter','Latex')
xlabel('t (s)')

Vmean=[mean(abs(diff(L)/dt)),mean(abs(diff(L2)/dt)),mean(abs(diff(L3)/dt)),mean(abs(diff(Lf)/dt))];
figure;
bar(Vmean)
set(gca, 'XTick',1:4)
set(gca,'XTickLabel',{['OU, \tau=',num2str(tau1),'s'],['OU, \tau=',num2str(tau2),'s'],['OU, \tau=',num2str(tau3),'s'],'LPOU, f_c=1Hz'})
ylabel('$\overline{|\dot{x}|}$','Interpreter','Latex')