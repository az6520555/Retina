clear
close all
Tot = 300;  
dt = 0.01; 
T = dt:dt:Tot;
D =4; %dynamical range
G=[20 10 6 3]
at = 60;
m = 10;
rand_seed=randn(1,length(T))
a2=[];ey=[];ey1=[];
path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\stimulus _data\20200419\';

for i=1:length(G)
    w = G(i)/(2*1.06);
    if i~=1
        Lold=L;
    end
    
    L = zeros(1,length(T));
    V = zeros(1,length(T));
    for t = 1:length(T)
        L(t+1) = L(t) + V(t)*dt;
        V(t+1) = (1-G(i)*dt)*V(t) - w^2*L(t)*dt + sqrt(dt*D)*rand_seed(t);
    end      
    
    if i==1
        Lold=L;
    end
    
    %% shift
    [val,lags] =xcorr(Lold,L,1000);
    [M,I]=max(val);
    figure(1);plot(lags,val)
    Lshift=L(-lags(I)+1:end);
    
    %% generate stimulus
    rate=20000;
    m=10;
    Lshift = m/5*Lshift/std(Lshift);
    temp = Lshift-mean(Lshift)+m; % X is the final stimuX
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
    figure(2);hold on;plot(t,ey);
    save([path,'HMM_G=',num2str(G(i)),'.mat'],'ey','a2','t')
end

