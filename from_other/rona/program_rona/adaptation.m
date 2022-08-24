% Analysis data that similar with "Coordinated dynamic encoding in the
% retina using opposing forms of plasticity", by Rona, July.2017
clear all
close all
clc

cd('E:\google_rona\20171020') ; % the folder of the files
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file) ; 

rr =[9,17,25,33,41,49,...
  2,10,18,26,34,42,50,58,...
  3,11,19,27,35,43,51,59,...
  4,12,20,28,36,44,52,60,...
  5,13,21,29,37,45,53,61,...
  6,14,22,30,38,46,54,62,...
  7,15,23,31,39,47,55,63,...
    16,24,32,40,48,56];
%%%%%%%%%%%%%%   user's setting  %%%%%%%%%%%%%%%%%%%%%%
BinningInterval=20e-3;
L1 = 2; %high contrast period
L2 = 15; %low contrast period
SamplingRate=20000;
trst=1;trnum=60; %sweep trail from trst to trend
f=10;p=1000;
roi=[1:60];

for m = 3
    clearvars -except all_file n_file m BinningInterval SamplingRate L1 L2 trst trnum f p roi rr
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-'; 
    TimeStamps2=TimeStamps(1:1:end); 
    if length(TimeStamps2)<=(trst+trnum-1)
        trnum = length(TimeStamps2)-trst+1;
    end
    trend = trst+trnum-1;
    DataTime = mean(diff(TimeStamps));  
    BinningTime = [ 0 : BinningInterval : DataTime];

    sweepend = trend;
    for sweepindex=trst:sweepend-1
        TimeStampsSweep=TimeStamps2(sweepindex:sweepindex+1); % forcus on which trails 
        cut_spikes = seperate_trials(Spikes,TimeStampsSweep);      
        for i = 1:60
            [n,xout] = hist(cut_spikes{i},BinningTime) ;
            BinningSpike(sweepindex,i,:) = n ;
        end
    end     
        %%%%%% a3 %%%%%%%%%
    y1=a_data(3,TimeStamps2(1)*20000:TimeStamps2(2)*20000);
    x1=1/20000:1/20000:length(y1)/20000;

        %%%%%%%%%%%%%%%%% raster plot %%%%%%%%%%%%%%%%%%%
    for nn = 1:length(roi)
        BinningSpike2 = sum(BinningSpike,1);
        SB = sum(BinningSpike2(1,roi(nn),:),2);
        SB1 = squeeze(SB)/size(BinningSpike,1)/BinningInterval;   

    %     figure(1)
    %     imagesc(BinningTime,[1:60],squeeze(BinningSpike2));
    %     title([name,'(sum over ',num2str(trend-trst),' trails) ']);  
    %     xlabel('t(s)');ylabel('channel');
    %     colorbar;
    %     saveas(gcf,[name,'_raster.jpg'],'jpg')
    %     saveas(gcf,[name,'_raster.fig'],'fig')    
    %%%%%%%%%%%%%%%%%%%%%%%%%  Plot histogram   %%%%%%%%%%%%%%%%%        
        figure(100);hold on;subplot(8,8,rr(roi(nn)));
%         figure(roi(nn));
        area(BinningTime(1:L1/BinningInterval),SB1(1:L1/BinningInterval),'FaceColor','b','EdgeColor','b');hold on; 
        area(BinningTime(L1/BinningInterval+1:(L1+10)/BinningInterval),SB1(L1/BinningInterval+1:(L1+10)/BinningInterval),'FaceColor','r','EdgeColor','r')
%         set(gca,'XTick',[])
    %%%%%%%%%%%%%%%%%%%%%%%%%  Adaptive index   %%%%%%%%%%%%%%%%% 
        re(roi(nn)) = sum(BinningSpike2(1,roi(nn),(L1+0.5)/BinningInterval:(L1+3.5)/BinningInterval));
        rl(roi(nn)) = sum(BinningSpike2(1,roi(nn),(L1+10)/BinningInterval:(L1+13)/BinningInterval));
        A(roi(nn)) = (re(roi(nn))-rl(roi(nn)))/(re(roi(nn))+rl(roi(nn)));
    end
    A(A==0)=NaN;figure;plot(A,'*')
end