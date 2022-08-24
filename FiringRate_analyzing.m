clear
close all
load('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200408\20200408_OU_cutoff=4_sort_unit1.mat')

channel=43;
Spikes{51}=[];
rate=20000;
t=1/rate:1/rate:length(a_data(1,:))/rate;
t=t(1:10:end); % downsample
a=a_data(1,:); a=a(1:10:end); % downsample
a(t<TimeStamps(1))=NaN;
a(t>TimeStamps(end))=NaN;
BinningInterval = 0.01;
BinningTime = [0 : BinningInterval : TimeStamps(end)];
BinningSpike = zeros(60,length(BinningTime));
for i = 1:60
    [n,xout] = hist(Spikes{i},BinningTime);
    BinningSpike(i,:) = n;
end
BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;

% sumBinningSpike=sum(BinningSpike,1);

yyaxis right
plot(BinningTime,BinningSpike(channel,:));hold on 
yyaxis left
plot(t,a)
