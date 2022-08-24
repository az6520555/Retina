%FFT of stimulus
clear 
% close all
cd(['G:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20190131\']);
all_file = dir('*.mat') ;
nnn=[7];
load(all_file(nnn).name)
SamplingRate=20000;
BinningInterval=0.01;

sti=a_data(3,:);
sti=sti(TimeStamps(1)*SamplingRate:TimeStamps(2)*SamplingRate);
sti_ave=sti-mean(sti);

[b,a] = butter(2,50/20000,'low'); % set butter filter
a_data2 = filter(b,a,a_data(3,:));
TriggerData = a_data2(TimeStamps(1)*SamplingRate:TimeStamps(length(TimeStamps))*SamplingRate);% figure;plot(isi);
inten = downsample(TriggerData,SamplingRate*BinningInterval);
inten2 = inten-mean(inten);
ttt=(1:length(inten))/100;

L=length(inten2);
y=fft(inten2);
f=(0:L/2)/L/BinningInterval;
P2=abs(y/L);
P1=P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(f,P1)

