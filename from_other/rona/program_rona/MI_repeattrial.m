% Calculate MI for repeat trail HMM, by Rona, Jun.2017
clear all
close all
cd('E:\google_rona\20161122\HMM_trail');
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
SamplingRate=20000;
cc=hsv(n_file);

for z =2
clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA ro
file = all_file(z).name ;
[pathstr, name, ext] = fileparts(file);
directory = [pathstr,'\'];
filename = [name,ext];
load([filename]);
name(name=='_')='-';
bin=10;  BinningInterval = bin*10^-3; 
backward=ceil(500/bin); forward=ceil(500/bin);
Nroi = 6;
Proi = 38;

% method = input('How to ultilize the spike trians? whole(w)/combine(c)/sum(s)')
method = 'sum';

if method == 'w'
    DataTime = (TimeStamps(end) - TimeStamps(1));
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(1)+DataTime];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end   
    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;
    s1 = BinningSpike(Nroi,:);
    s2 = BinningSpike(Proi,:);

elseif method == 'c'
    BinningSpike2 = []; adata = [];
    for j = 1:20
        DataTime = (TimeStamps(2) - TimeStamps(1))/2;
        BinningTime = [TimeStamps(j) : BinningInterval : TimeStamps(j)+DataTime];
        BinningSpike = zeros(60,length(BinningTime));
        for i = 1:60
            [n,xout] = hist(Spikes{i},BinningTime) ;
            BinningSpike(i,:) = n ;
        end
        BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;
        BinningSpike2 = [BinningSpike2 BinningSpike];
    end
    s1 = BinningSpike2(Nroi,:);
    s2 = BinningSpike2(Proi,:);
 
else
    DataTime = (TimeStamps(2) - TimeStamps(1))/2;
    BinningTime = [0 : BinningInterval :  DataTime];
    sweepend=length(TimeStamps);
    set(gcf,'position',[150,30,1024,900])
    h = subplot(sweepend,1,1);
    for sweepindex=1:sweepend-1
        TimeStampsSweep=TimeStamps(sweepindex:sweepindex+1); % forcus on which trails 
        cut_spikes = seperate_trials(Spikes,TimeStampsSweep);      
        for i = 1:60
            [n,xout] = hist(cut_spikes{i},BinningTime) ;
            BinningSpike(sweepindex,i,:) = n ;
        end
        BinningSpike(:,:,1) = 0;BinningSpike(:,:,end) = 0;
        subplot(sweepend,1,sweepindex);hold on
        plot(BinningTime,squeeze(sum(BinningSpike(sweepindex,Proi,:),2)),'r');
        plot(BinningTime,squeeze(sum(BinningSpike(sweepindex,Nroi,:),2)),'b');
    end
    samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
    BinningSpike2 = squeeze(sum(BinningSpike(:,:,:),1));
    s1 = BinningSpike2(Nroi,:);
    s2 = BinningSpike2(Proi,:);
end

%% Stimulation
isi = a_data(3,TimeStamps(1)*20000:(TimeStamps(1)+DataTime)*20000);
isi2 = [];
states = 8;
X = isi;
nX = sort(X);
abin = length(nX)/states;
intervals = [nX(1:abin:end) inf]; 
temp=0;
for jj = 1:BinningInterval*SamplingRate:length(X)
    temp=temp+1;
    isi2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
end
if method == 'c'
    isi2 = repmat(isi2,1,20);
else
    isi2 = isi2;
end

%% delete repeat spikes
%  for i=12:length(s2)
%      if s2(i) ~= 0 
%          s1(i-110/bin:i-90/bin)=0;
%      else
%      end
%  end
 
 %% faster s2 spikes
%  for i=12:length(s2)
%      if s2(i) ~= 0 
%          s2(i-100/bin) = s2(i);
%          s2(i) = 0;
%      else
%      end
%  end

%% mutual information
 Neurons = s1; %N
 isi2 = s2; %P
    %% shuffle
%         r=randperm(length(Neurons));
%         for j=1:length(r)            
%             sNeurons(j)=Neurons(r(j));
%         end
%         Neurons=sNeurons;

dat=[];informationp=[];temp=backward+2;
for i=1:backward+1 %past(t<0)
    x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
    y = isi2(forward+1:length(isi2)-backward)';
    dat = [x,y];
    
    [N,C]=hist3(dat,[max(Neurons)+1,8]); %20:dividing firing rate  6:# of stim
    px=sum(N,1)/sum(sum(N)); % x:stim
    py=sum(N,2)/sum(sum(N)); % y:word
    pxy=N/sum(sum(N));
    temp2=[];
    for j=1:length(px)
        for k=1:length(py)
          temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
        end
    end
    temp=temp-1;
    informationp(temp)=nansum(temp2(:));
    c=corrcoef(x,y);
    corrp(temp)=c(2,1);
end  

dat=[];informationf=[];temp=0;sdat=[];
for i=1:forward
    x =Neurons(forward+1-i:length(Neurons)-backward-i)';
    y = isi2(forward+1:length(isi2)-backward)';
    dat=[x,y];

    [N,C]=hist3(dat,[max(Neurons)+1,8]); %20:dividing firing rate  6:# of stim
    px=sum(N,1)/sum(sum(N)); % x:stim
    py=sum(N,2)/sum(sum(N)); % y:word
    pxy=N/sum(sum(N));
    temp2=[];
    for j=1:length(px)
        for k=1:length(py)
            temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
        end
    end
    temp=temp+1;
    informationf(temp)=nansum(temp2(:)); 
    c=corrcoef(x,y);
    corrf(temp)=c(2,1);
end

information = [informationp informationf];
fr(z)=mean(Neurons)/BinningInterval; %firing rate(Hz)

t=[-backward*bin:bin:forward*bin];  
figure(110);hold on;plot(t,information,'LineWidth',2,'DisplayName',num2str(n));%'color',cc(z,:),
% legendInfo{z} = [num2str(z)];legend(legendInfo);
% legendInfo{nn} = ['#' num2str(n)];legend(legendInfo);
set(gca,'FontSize',16)
xlabel('\deltat (ms)');ylabel('MI (bits)');
end
