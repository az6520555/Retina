%%% calculate correlation between firing among different neurons, Rona, March.2016 %%%
clear all
% close all
clc

cd('\\192.168.0.100\Experiment\Retina\Rona\Exp\20160229') ; % the folder of the files
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file) ; 

%%%%%%%%%%%%%%   user's setting  %%%%%%%%%%%%%%%%%%%%%%
BinningInterval=0.005;
SamplingRate=20000;

for m = 1
    clearvars -except x y all_file n_file m BinningInterval SamplingRate chst1 chend1 chst2 chend2 trst trnum
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';

%%%%%%%%%%%%%%%%%%%%%%%%%%  Binning  %%%%%%%%%%%%%%%%%%%%%%%%
    DataTime = (TimeStamps(3) - TimeStamps(2));
    BinningTime = [ 0 : BinningInterval : DataTime];
   
    CC=[];
    sweepend=2;
    
    for k =1:10
        corst = 10; %calculate from 10s withing a trail
        corend = 0.5; %caculating period (unit:s) 
        t1 = corst/BinningInterval+1+(k-1)*corend/BinningInterval;
        t2 = t1 + corend/BinningInterval;
       
        for sweepindex = 1:length(TimeStamps)-1
            TimeStampsSweep=TimeStamps(sweepindex:sweepindex+1); % forcus on which trails 
            cut_spikes = seperate_trials(Spikes,TimeStampsSweep);      
            for i = 1:60
                [n,xout] = hist(cut_spikes{i},BinningTime) ;
                BinningSpike(i,:) = n ;
            end
            corrBinningSpike = BinningSpike(:,t1:t2);
            corrBinningSpike(find(BinningSpike(:,t1:t2)>1)) = 1;

            for i = 1:size(corrBinningSpike,1)
                for j = 1:size(corrBinningSpike,1)
                    if corrBinningSpike(i,:)==0 & corrBinningSpike(j,:)~=0
                        A=0;
                    else
                        A = corrcoef(corrBinningSpike(i,:),corrBinningSpike(j,:));
                    end
                    B(i,j) = A(1,2);  %correlation of all neurons                   
                end
            end
            CC(sweepindex,:,:) = B;
        end
        M(k,:,:) = nanmean(CC,1);
        C = reshape(M(k,:,:),60,60);
        
        figure(k);imagesc(C);
        title([num2str(t1*BinningInterval),'s to ',num2str(t2*BinningInterval),'s']);
        colorbar;
        
        for i = 1:60
            fr(k,i) = sum(BinningSpike(i,t1:t2)) ;
        end
        colorbase={[1,0,0],[0,1,0],[0,0,1],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]};
        t = (corend:corend:k*corend);
        figure(1);hold on
        plot(t,fr(:,39),'color',colorbase{m});
%         legend(name);hold on
        title(['Firing Rate']);xlabel('time after stimuli(s)');ylabel(['FiringRate per channel per ',num2str(corend),'s']);
        
        figure(2);plot(t,M(:,39,45),'color',colorbase{m});hold on;
        title(['Correlation']);xlabel('time after stimuli(s)');ylabel('Correlation');
  
    end

end