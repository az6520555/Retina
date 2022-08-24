% Calculate fft for a spike train, by Rona 
clear all
close all
load('E:\google_rona\20171128\spon_20171128.mat');
BinningInterval=0.001;
SamplingRate=20000;
TimeStamps=[10,120]; % for spontaneous
TimeStamps2=TimeStamps(1:1:length(TimeStamps)); 
cut_spikes = Spikes;    

BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(2)];


%%%%%%%%%%%%%%%%%%%%%%%%%  Plot Different Trials   %%%%%%%%%%%%%%%%% 
figure(2);
for i = 1:60
    [n,xout] = hist(Spikes{i},BinningTime) ;
    BinningSpike(i,:) = n ;
    BinningSpike(i,1)=0; BinningSpike(i,end)=0;
end
BinningSpike2 = sum(BinningSpike,1);
plot(BinningTime,BinningSpike2);

%%%%%%%%%%%%%%%%% raster plot %%%%%%%%%%%%%%%%%%%
figure(1)
imagesc(BinningTime,[1:60],BinningSpike);
colorbar;

%% fft1
X = BinningSpike2;
NFFT = 2^nextpow2(size(X,2));
% Y = fft(smooth(BinningSpike,kkk,'loess'),NFFT) / size(smooth(BinningSpike,kkk,'loess'),1);  %fft(BinningSpike(i,:),NFFT)/size(BinningSpike,2);
Y = fft(X,NFFT)/size(X,2);
Fs=1/BinningInterval;
f = Fs/2*linspace(0,1,NFFT/2+1);
[c, index] = min(abs(f-2));
figure; plot(f,2*abs(Y(1:NFFT/2+1)),'K');

%% fft2
% X=BinningSpike2;
Y = fft(X);
L=length(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs=1/BinningInterval;
f = Fs*(0:(L/2))/L;
figure;plot(f,P1)
