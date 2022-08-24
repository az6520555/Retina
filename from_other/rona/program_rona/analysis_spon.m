clear all
close all
clc

cd('\\192.168.0.100\experiment\Neuron\YiKo\checker') ; % the folder of the files
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file) ; 

%%%%%%%%%%%%%%   user's setting  %%%%%%%%%%%%%%%%%%%%%%
BinningInterval=0.001;
SamplingRate=20000;
chst1=1;chend1=60;%focus on channel chst1= to chend1
chst2=38;chend2=45;
trst=1;trnum=10; %sweep trail from trst to trend
f=10;p=20;

for m =3
    clearvars -except all_file n_file m BinningInterval SamplingRate chst1 chend1 chst2 chend2 trst trnum f p ch
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%  Binning  %%%%%%%%%%%%%%%%%%%%%%%%
    TimeStamps=[0:300:600]; % for spontaneous
    TimeStamps2=TimeStamps(1:1:length(TimeStamps)); 
    if length(TimeStamps2)<=(trst+trnum-1)
        trnum = length(TimeStamps2)-trst+1;
    end
    trend=trst+trnum-1;
 
%     DataTime = (TimeStamps2(2) - TimeStamps2(1));
DataTime = 600;
    cut_spikes = seperate_trials(Spikes,TimeStamps2(trst:trend));    

    BinningTime = [ 0 : BinningInterval : DataTime];
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%  Plot Different Trials   %%%%%%%%%%%%%%%%% 
%     if length(TimeStamps2)<=18
%         sweepend=length(TimeStamps2);
%     else
%         sweepend=24;
%     end
    sweepend=trend;
    figure(2);
    set(gcf,'position',[150,30,1024,900])
    h = subplot(sweepend+1,1,1);
    
    for sweepindex=1:sweepend-1
        TimeStampsSweep=TimeStamps2(sweepindex:sweepindex+1); % forcus on which trails 
        cut_spikes = seperate_trials(Spikes,TimeStampsSweep);      
        for i = 1:60
            [n,xout] = hist(cut_spikes{i},BinningTime) ;
            BinningSpike(sweepindex,i,:) = n ;
        end

        subplot(sweepend+1,1,sweepindex);
%         plot(BinningTime,squeeze(sum(BinningSpike(sweepindex,ch,:),2)));
        plot(BinningTime,squeeze(sum(BinningSpike(sweepindex,chst1:chend1,:),2)));
    end
    
%     samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
    set(get(h,'title'),'string',[name,'  ch',num2str(chst1),' to ',num2str(chend1)]);

%     saveas(gcf,[name,'_trials.jpg'],'jpg');
%     saveas(gcf,[name,'_trials.fig'],'fig');  

 %%%%%%%%%%%%%%%%% raster plot %%%%%%%%%%%%%%%%%%%
     BinningSpike2 = sum(BinningSpike(2,:,:),1);
%     BinningSpike2 = sum(BinningSpike(trst:trend-1,:,:),1);
%     SB=sum(BinningSpike2(1,ch,:),2);
    SB=sum(BinningSpike2(1,chst1:chend1,:),2);
    SB1=squeeze(SB);   
 
    figure(1)
    imagesc(BinningTime,[1:60],squeeze(BinningSpike2));
     title([name,'  (the 2nd trail) ']);  
%     title([name,'(sum over ',sprintf('%.0f',length(TimeStamps2)),' trails) ']);  
    xlabel('t(s)');ylabel('channel');
    colorbar;
%     saveas(gcf,[name,'_raster.jpg'],'jpg')
%     saveas(gcf,[name,'_raster.fig'],'fig')    
    
%%%%%%%%%%%%%  get peaks  %%%%%%%%%%%%%%%%%%%%%   
% %     SB1=sum(BinningSpike(chst1:chend1,:),1)/trnum;  
%     [spks,slocs]=findpeaks(SB1,'minpeakdistance',floor(DataTime/BinningInterval/f),'MINPEAKHEIGHT',p);    
% 
%     if isempty(slocs)
%         slocstime=NaN;
%     end
%     slocstime=(slocs-1)*BinningInterval; 
     
%%%%%%%%%%%%%%%%%%%%%%%%%  Plot histogram   %%%%%%%%%%%%%%%%%        
%     figure(3)
%     subplot(6,1,2:5); 
%     plot(BinningTime,SB1,slocstime,spks,'r*'); 
%     title([name,'   Ch',num2str(chst1),' to ',num2str(chend1),sprintf('\n'),' PeakTime(s)=',sprintf('%8.3f',slocstime),sprintf('\n')]);
% %     ylabel('firing rate per 5ms');
%     subplot(6,1,6);
%     plot(x1,y1,'r');
%     plot(Time(1:length(TimeStamps(1)*20000:TimeStamps(2)*20000)),1800+Filter_a_data(TimeStamps(1)*20000:TimeStamps(2)*20000),'r');
%     t2=TimeStamps2(2);t1=TimeStamps2(1)
%     time=[1/SamplingRate:1/SamplingRate:length(a_data(round(t1*20000):round(t2*20000)))/20000];
   
    
%     plot(time,1/700*a_data(round(t1*20000):round(t2*20000)),'r');
%     plot(BinningTime,SB1,slocstime,spks,'r*');
%     axis([0 DataTime 0 max(SB1)]);
%     xlabel('t(s)');
   
%     saveas(gcf,[name,'_hist.jpg'],'jpg')
%     saveas(gcf,[name,'_hist.fig'],'fig')
               

%     save(filename,'x','y','chst1','chend1','slocstime', '-append'); % save data into original dat

    
%     catch
%         [msgstr,msgerr] = lasterr;
%         disp([msgstr,msgerr])
%     end
       
%%%%%%%%%  focus on other channels   %%%%%%%%%%%%%%    
%     clearvars spks slocs slocstime theta 
%     SB2=sum(BinningSpike(c:d,:))/length(TimeStamps);    
%     [spks,slocs]=findpeaks(SB2,'minpeakdistance',floor(length(BinningTime)/4),'MINPEAKHEIGHT',0.2);   
%     slocstime=(slocs-1)*BinningInterval;    
%     theta=slocs/length(BinningTime)*2; %firing phase(unit:pi)
%      wf=(y(slocs(1)+1)-y(slocs(1)-1))/(2*BinningInterval)/((y(slocs(1)+1)+y(slocs(1)-1))/2);  %%weber's fraction
%     figure
%     plot(BinningTime,SB2,slocstime,spks,'r*');
%     title(['T=', sprintf('%.3f',DataTime),'s','    Ch',num2str(c),' to ',num2str(d),sprintf('\n'),' PeakTime(s)=',sprintf('%8.3f',slocstime),sprintf('\n'),...
%         'PeakPhase(pi)=',sprintf('%8.3f',theta),sprintf('\n'),'Weber Fraction=',sprintf('%.3f',wf)]);
%     xlabel('t(s)');ylabel('spike# in 60 channels within 5ms');
%     saveas(gcf,['T=', num2str(DataTime),'s',' Ch',num2str(a),'to',num2str(b),'_hist.jpg'],'jpg')
%     saveas(gcf,['T=', num2str(DataTime),'s',' Ch',num2str(a),'to',num2str(b),'_hist.fig'],'fig')  
        
end