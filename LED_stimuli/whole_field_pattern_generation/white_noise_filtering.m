clear
close all
Tot = 300;  
dt = 0.01; 
T = dt:dt:Tot;
tau=1;
D =4; %dynamical range
at = 60;
m = 10;
a2=[];ey=[];ey1=[];
L = zeros(1,length(T));
L=randn(1,length(L));
path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\stimulus _data\20200506\';

cutoff=[7,4,2,1,0.5];
corrtimef=[];
for afc=1:length(cutoff)
    [b,a]=butter(2,cutoff(afc)*dt,'low')
    Lf=filter(b,a,L);

%     [val,lags] =xcorr(L,Lf,1000);
%     [M,I]=max(val);
%     figure(1);plot(lags,val)
%     Lfshift=Lf(-lags(I)+1:end);

    figure(2);hold on
    plot(T,Lf)
   

    %%  auto correlation
    aaa=[];
    lags=[];
    [aaa,lags]=autocorr(L,'NumLags',1000);
    [M,ind]=min(abs(aaa-0.5));
    corrtime=ind*dt;

    aaa=[];
    lags=[];
    [aaa,lags]=autocorr(Lf,'NumLags',1000);
    [M,ind]=min(abs(aaa-0.5));
    corrtimef=[corrtimef ind*dt];

    %% fft
    Fs = 1/dt
    Ts =  dt            % Sampling period       
    Len = length(L);             % Length of signal
    t = (0:Len-1)*Ts;        % Time vector

%     Y = fft(L);
%     P2 = abs(Y/Len);
%     P1 = P2(1:Len/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = Fs*(0:(Len/2))/Len;
%     figure(3);hold on;plot(f,smooth(P1,100),'linewidth',2) 
%     xlim([0 5])

    Yf = fft(Lf);
    P2 = abs(Yf/Len);
    P1 = P2(1:Len/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(Len/2))/Len;
    figure(3);hold on;plot(f,P1,'linewidth',2) 
    xlim([0 5])
    
    %% generate stimulus
    rate=20000;
    m=10;
    Lf = m/5*Lf/std(Lf);
    temp = Lf-mean(Lf)+m; % X is the final stimuX
    X = abs(temp);
    std(X)/m
    ey1=zeros(1,length(X));
    temp = X(2:1:length(X));%temp = X(2:dtau/dt:Tot/dt);
    temp2=repmat(temp,rate*dt,1);
    ey1=temp2(:)';

    ey0=m*ones(1,at*rate);%REST
    ey=[ey0 ey1];
    a2=zeros(1,length(ey));
    a2(at*rate:(at+1)*rate)=1;
    a2((Tot+at-1)*rate:(Tot+at)*rate)=1;
    t=[1/rate:1/rate:length(ey)/rate];
    figure(4);hold on;plot(t,ey);
%     save([path,'WhiteNoise_cutoff=',num2str(cutoff(afc)),'.mat'],'ey','a2','t')
end
figure(2);plot(T,L);
L = m/5*L/std(L);
temp = L-mean(L)+m; % X is the final stimuX
X = abs(temp);
std(X)/m
ey1=zeros(1,length(X));
temp = X(2:1:length(X));%temp = X(2:dtau/dt:Tot/dt);
temp2=repmat(temp,rate*dt,1);
ey1=temp2(:)';

ey0=m*ones(1,at*rate);%REST
ey=[ey0 ey1];
a2=zeros(1,length(ey));
a2(at*rate:(at+1)*rate)=1;
a2((Tot+at-1)*rate:(Tot+at)*rate)=1;
t=[1/rate:1/rate:length(ey)/rate];
figure(4);hold on;plot(t,ey);
% save([path,'WhiteNoise_.mat'],'ey','a2','t')