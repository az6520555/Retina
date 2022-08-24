% Calculate information amount with various number of neurons, revised from
% Kevin's code by Rona

clear all
close all

cd('G:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20190329');
all_file = dir('*.mat'); % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
for z = [10]

SamplingRate=20000;
cc=hsv(n_file);
file = all_file(z).name ;
[pathstr, name, ext] = fileparts(file);
directory = [pathstr,'\'];
filename = [name,ext];
load([filename]);
SamplingRate=20000;

%% Stimuli process
isi = a_data(3,TimeStamps(1)*20000: TimeStamps(length(TimeStamps))*20000);
DataTime = max(TimeStamps); Start = min(TimeStamps);
bin = 10;  BinningInterval = bin*10^-3;  %ms
shift = [500 500];  %backward and forward shiftin ms
resevoir = [43 44 46];%3,6,8,10,15,18,30,56

for rrr = 1:1:length(resevoir)
    c = combnk(1:length(resevoir),rrr);
    
for rr = 1:size(c,1)
    roi = resevoir(c(rr,:));

%% Spike procss
BinningTime = [Start : BinningInterval : DataTime];
BinningSpike = zeros(length(roi),length(BinningTime));
SpikeLate = zeros(size(BinningSpike));
for i = 1:length(roi)
    temp = Spikes{roi(i)};
    temp = temp(find(Start<temp));
    [n,xout] = hist(temp,BinningTime) ;
    BinningSpike(i,:) = n;
end

isi2=[];
states=8;
X=isi;
nX = sort(X);
abin = length(nX)/states;
intervals = [nX(1:abin:end) inf]; 
temp=0;
for jj = 1:BinningInterval*SamplingRate:length(X)
    temp=temp+1;
    isi2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
    inten(temp)= X(jj);
end

%% Mutual information
Neurons = sum(BinningSpike,1);

    %%%randomize
%     randomness = randperm(length(Neurons));
%     Neurons = Neurons(randomness);
%     isi2 = isi2(randomness);

backward=round(shift(1)/bin); forward=round(shift(2)/bin);
dat=[];informationp=[];temp=backward+2; 
for i=1:backward+1
    dat=[Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))',Neurons(forward+1:length(isi2)-backward)'];
    [N,C]=hist3(dat,[floor(max(Neurons)+1),max(isi2)]); %20:dividing firing rate  6:# of stim
    px=sum(N,1)/sum(sum(N)); % x:stim
    py=sum(N,2)/sum(sum(N)); % y:word
    pxy=N/sum(sum(N));
    temp2=[];
    nor = sum(Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1)));% / length(Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1)));
    for j=1:length(px)
        for k=1:length(py)
          temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);%  / nor;  %
        end
    end
    temp=temp-1;
    informationp(temp)=nansum(temp2(:)); 
end   

dat=[];informationf=[];temp=0;sdat=[];
for i=1:forward
    dat=[Neurons(forward+1-i:length(Neurons)-backward-i)',Neurons(forward+1:length(isi2)-backward)'];
    [N,C]=hist3(dat,[floor(max(Neurons)+1),max(isi2)]); %20:dividing firing rate  6:# of stim   %max(Neurons)+1,max(isi2)
    px=sum(N,1)/sum(sum(N)); % x:stim
    py=sum(N,2)/sum(sum(N)); % y:word
    pxy=N/sum(sum(N));
    temp2=[];
    nor = sum(Neurons(forward+1-i:length(Neurons)-backward-i));% / length(Neurons(forward+1-i:length(Neurons)-backward-i));
    for j=1:length(px)
        for k=1:length(py)
            temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);%  / nor;  %
        end
    end
    temp=temp+1;
    informationf(temp)=nansum(temp2(:)); 
end
information=[informationp informationf]/ BinningInterval;

t=[-backward*bin:bin:forward*bin];
% figure(1);hold on; 
% plot(t,information,'color',[.47 .67 .19],'LineWidth',2.5); set(gca,'FontSize',16)
% % title(['Neuron',num2str(n),],'FontSize',10);
% xlabel('\delta t(ms)','FontSize',20);ylabel('mutual information','FontSize',20);
% figure;plot(isi2);hold on;plot(Neurons,'*')

    a=1;b=length(information);
    [spks,slocs]=findpeaks(information(a:b),'minpeakdistance',b-a-1);
    slocstime=(slocs-1)*bin-backward*bin+bin*(a-1);
%     w(nn,z)=widths(1)*bin;
    statispks(rr,rrr)=spks;
    statistime(rr,rrr)=slocstime;  
end
end
    statispks(statispks==0)=NaN;
    figure;plot(nanmean(statispks))
end

