clear all
close all

path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\temporal_pattern\20181220\';
cd(path)
load('20181220_OSR_ChangeFrequency_200to210.mat')
SamplingRate=20000;
a_data3=a_data(3,:);
t=1/SamplingRate:1/SamplingRate:length(a_data(3,:))/SamplingRate;
indexofTimeStamps=1:40:length(TimeStamps);
BinningInterval=0.01;
BinningTime= [0 : BinningInterval : TimeStamps(indexofTimeStamps(2))-TimeStamps(indexofTimeStamps(1))];


stimulus=cell(1,length(indexofTimeStamps)-1)
spikessum=cell(1,60);
for i=1:length(indexofTimeStamps)-1
    timevector=[TimeStamps(indexofTimeStamps(i)) TimeStamps(indexofTimeStamps(i+1))];
    for j=1:60
        spikeST=[];
        spikeST=Spikes{j}(Spikes{j}>timevector(1) & Spikes{j}<timevector(2))-timevector(1);
        spikessum{j}=[spikessum{j} spikeST];
    end
end
BinningSpikes=cell(1,60);
for k=1:60
    [n,xout] = hist(spikessum{k},BinningTime);
    BinningSpikes{k}=n;
end

stimulus=a_data3(t>TimeStamps(indexofTimeStamps(1)) & t<TimeStamps(indexofTimeStamps(2)));
tsti=1/SamplingRate:1/SamplingRate:length(stimulus)/SamplingRate;
subplot(2,1,1)
plot(BinningTime,BinningSpikes{1})
subplot(2,1,2)
plot(tsti,stimulus)
samexaxis('abc','xmt','on','ytac','join','yld',1);

for i=1:60
    if isempty(spikessum{i})==1
        spikessum{i}=0;
    end
end
figure(87);
subplot(10,1,[1:9])
LineFormat.Color = [0.1 0.1 0.1];
plotSpikeRaster(spikessum,'PlotType','vertline','LineFormat',LineFormat)
subplot(10,1,10)
plot(tsti,stimulus)
samexaxis('abc','xmt','on','ytac','join','yld',1);