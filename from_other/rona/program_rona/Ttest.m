% Calculated P value between N and P cells, by Rona, Dec.2017
clear all
% close all
clc
BinningInterval=0.005;
SamplingRate=20000;
load('E:\google_rona\20170504\onoff_before_20170504.mat') ; % the folder of the files

%%   user's setting  %%%%%%%%%%%%%%%%%%%%%%
chst1=1;chend1=60;%focus on channel chst1= to chend1
trst=1;trnum=20; %sweep trail from trst to trend
f=10;p=100;

Proi=[1,6,12,18,23,24,25,26];
Nroi=[2,3,10,11,15,16,19,31,32,46,51,57];

%%  Binning  %%%%%%%%%%%%%%%%%%%%%%%%
    TimeStamps2=TimeStamps(1:4:length(TimeStamps)); 
    if length(TimeStamps2)<=(trst+trnum-1)
        trnum = length(TimeStamps2)-trst+1;
    end
    trend=trst+trnum-1;
 
    DataTime = (TimeStamps2(2) - TimeStamps2(1));

    cut_spikes = seperate_trials(Spikes,TimeStamps2(trst:trend));    

    BinningTime = [ 0 : BinningInterval : DataTime];
  
    sweepend=trend;
    
    for sweepindex=1:sweepend-1
        TimeStampsSweep=TimeStamps2(sweepindex:sweepindex+1); % forcus on which trails 
        cut_spikes = seperate_trials(Spikes,TimeStampsSweep);      
        for i = 1:60
            [n,xout] = hist(cut_spikes{i},BinningTime) ;
            BinningSpike(sweepindex,i,:) = n ;
        end
    end     
BinningSpike2 = squeeze(sum(BinningSpike(trst:trend-1,:,:),1));
   
%% sustained
norm = 1/BinningInterval;
sus = mean([sum(BinningSpike2(:,2*norm:4*norm),2),sum(BinningSpike2(:,6*norm:8*norm),2),sum(BinningSpike2(:,10*norm:12*norm),2)],2)/2;

%% MI
% clearvars -except sus SamplingRate Proi Nroi
load('E:\google_rona\20170504\HMM\HMM_G=10_20170504.mat');
roi = [1:60];
bin=10;  BinningInterval = bin*10^-3; 
backward=ceil(500/bin); forward=ceil(500/bin);

    %% a_data as isi  %   a_data2=(a_data-32768)*0.1042;
    [b,a] = butter(2,50/20000,'low'); % set butter filter
    a_data2 = filter(b,a,a_data(3,:));
    isi = a_data2(TimeStamps(1)*20000:TimeStamps(length(TimeStamps))*20000);% figure;plot(isi);
 %% Spike process
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end
    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;
    %% state of light intensity %%% 
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

%% Mutual Information
    for nn=1:length(roi)
        n=roi(nn);
        Neurons = BinningSpike(n,:); % Neurons = sum(BinningSpike(:,:));
         %% shuffle
%         r=randperm(length(Neurons));
%         for j=1:length(r)            
%             sNeurons(j)=Neurons(r(j));
%         end
%         Neurons=sNeurons;

    dat=[];informationp=[];temp=backward+2;
    for i=1:backward+1 %past(t<0)

        if length(Neurons)==length(isi2)
            x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
        else
            x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1)-1)';
        end    
        y=isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        [N,C]=hist3(dat{i},[max(Neurons)+1,8]); %20:dividing firing rate  6:# of stim
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
        if length(Neurons)==length(isi2)
            x =Neurons(forward+1-i:length(Neurons)-backward-i)';
        else
            x = Neurons(forward+1-i:length(Neurons)-backward-i-1)';
        end     
        y = isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        [N,C]=hist3(dat{i},[max(Neurons)+1,8]); %20:dividing firing rate  6:# of stim
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

    information=[informationp informationf]/BinningInterval;
    corr=[corrp corrf];
%     fr(z)=mean(Neurons)/BinningInterval; %firing rate(Hz)
    [pks(1,nn),plocs(1,nn)]=max(information);
    t=[-backward*bin:bin:forward*bin];  
  dtp = t(plocs)';
    end
figure;plot(dtp(Proi),sus(Proi),'*')
hold on;plot(dtp(Nroi),sus(Nroi),'r*')
figure;bar([mean(sus(Nroi)),mean(sus(Proi))]);hold on
set(gca, 'FontSize',12,'XTick',[1 2],'XTickLabel',{'N cells','P cells'});
errbar = [std(sus(Proi)) std(sus(Proi))];
yd = [mean(sus(Nroi)),mean(sus(Proi))];
hold on;errorbar([1:2],  yd,  errbar, '.k', 'LineWidth',2);hold off