clear
close all
itau=[0.5,0.8,1.5,3]
for j =1
    Tot = 1000;  
    dt = 0.01; 
    T = dt:dt:Tot;
    tau=itau(j);
    D =4; %dynamical range
    at = 120;
    m = 10;
    a2=[];ey=[];ey1=[];
    L = zeros(1,length(T));
    for i = 1:length(T)-1
        L(i+1) = L(i) + (-L(i)/tau + randn*sqrt(D/dt))*dt;
    end  
    % L = m*L/(std(L)*5);
    % X = L-mean(L)+m; % X is the final stimuX
    % X = abs(X);

%     [b,a]=butter(2,0.3*dt,'low')
%     Lsm=filter(b,a,L);
%     R=Lsm-L;
    
    Lsm=smooth(L,21);
    R=Lsm'-L;
    
    figure(1);hold on
    plot(T,L)
    plot(T,Lsm)
%     plot(T,R)

    % for i=1:length(X)
    %     ey1(rate*dt*(i-1)+1:rate*dt*i)=X(i);
    % end
    % std(X)/m
    % ey0=m*ones(1,at*rate);%REST
    % ey=[ey0 ey1];
    % a2=zeros(1,length(ey));
    % a2(at*rate:(at+1)*rate)=1;
    % a2((Tot+at-1)*rate:(Tot+at)*rate)=1;
    % t=[1/rate:1/rate:length(ey)/rate];
    % figure(3);plot(t,ey);hold on 

    aaa=[];
    lags=[];
    [aaa,lags]=autocorr(L,'NumLags',1000);
    [M,ind]=min(abs(aaa-0.5));
    corrtime=ind*dt;

    [aaa,lags]=autocorr(Lsm,'NumLags',1000);
    [M,ind]=min(abs(aaa-0.5));
    corrtime_sm=ind*dt;

    %% fourier

    Fs = 1/dt
    Ts =  dt            % Sampling period       
    Len = length(L);             % Length of signal
    t = (0:Len-1)*Ts;        % Time vector

    Y = fft(L);
    P2 = abs(Y/Len);
    P1 = P2(1:Len/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(Len/2))/Len;
    figure(2);hold on;plot(f,P1) 
    xlim([0 5])

    Ysm = fft(Lsm);
    P2 = abs(Ysm/Len);
    P1 = P2(1:Len/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(Len/2))/Len;
    figure(2);hold on;plot(f,P1) 
    xlim([0 5])

    Ysm = fft(R);
    P2 = abs(Ysm/Len);
    P1 = P2(1:Len/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(Len/2))/Len;
    figure(4);hold on;plot(f,P1) 
end
