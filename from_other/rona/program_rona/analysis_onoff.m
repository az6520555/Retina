
% Analysis onoff experiment, by Rona
clear all
close all
clc
path_sort=['I:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\sorted data\20190117\deadtime_3000\'];
try
cd('\\nasd32ad0\experiment\Public\ToChou\Leo0807exp\new') ; % the folder of the files
catch
cd('') ; % the folder of the files    
end
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file) ; 
%%%%%%%%%%%%%%   user's setting  %%%%%%%%%%%%%%%%%%%%%%
BinningInterval=0.01;
SamplingRate=20000;
roi=[1:60]; %region of interest
trst=1;trnum=7; %sweep trail from trst to trend
f=10;p=100;

for m = [6]
    clearvars -except all_file n_file m BinningInterval SamplingRate chst1 chend1 chst2 chend2 trst trnum f p ch roi path_sort
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);

    Spikes=Spikes(1,:);
    
%     name(name=='_')='-';
 
%     load(['D:\rona\Google Drive\20160427\diode.mat']);
%     load(['C:\Users\HydroMatlab\Google ¶³ºÝµwºÐ\20160429\diode\','diode_',filename]);
%%%%%%%%%%%%%%%%%%%%%%  TimeStamps  %%%%%%%%%%%%%%%%%%%
    a_data2 = a_data(2,:);
%     a_data2 = a_data2 - a_data2(1);
%     a_data = diff(a_data);
    [~,locs]=findpeaks(diff(a_data2),'MINPEAKHEIGHT',5*std(diff(a_data2)));
     TimeStamps = (locs)/SamplingRate;
%% Spike process
%     xls=xlsread([name,'.xls']);
%     ss = [29,18,39,53,4,31,59,5,23,38,13,15,...
%         42,37,21,17,55,35,28,9,47,54,10,34,...
%         60,11,32,43,12,25,57,20,16,56,38,22,...
%         8,48,46,30,2,40,44,3,33,52,19,24,...
%         51,6,26,50,14,7,49,36,27,1,41,45];
%     Spikes=cell(1,60);
%     for j = 1:max(xls(:,1))
%             temp = 0;
%         for i = 1:length(xls)
%             if xls(i,1) == j
%                 temp = temp+1;
%                 Spikes{ss(j)}(1,temp) = xls(i,2);
%             end
%         end
%     end    
% %================== new sorting sequence ========================
%     ss=[29,30,28,27,22,21,14,20,...
%         13,6,12,5,19,11,4,3,10,...
%         18,2,9,1,8,17,7,16,15,26,...
%         25,23,24,32,31,33,34,39,...
%         40,47,41,48,55,49,56,42,...
%         50,57,58,51,43,59,52,60,...
%         53,44,54,45,46,35,36,38,37];
% 
% %================= for sorted data in .mat files ==========================
%     Spikes = cell(1,60);
%     load([path_sort,filename(1:end-4),'_sort.mat']);
%     temp_spikes={};
%     for h=1:60
%         if h<11
%             temp_spikes{ss(h)} = eval(['adc00',int2str(h-1)]);
%         else
%             temp_spikes{ss(h)} = eval(['adc0',int2str(h-1)]);
%         end
%     end
%     for i=1:60
%         if isempty(temp_spikes{i})==1
%             continue
%         end
%         for j=1:length(temp_spikes{i}(:,1))
%             if temp_spikes{i}(j,3)==1 % this determine which unit we choose
%                 Spikes{i}=[Spikes{i} temp_spikes{i}(j,1)];
%             end
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%  Binning  %%%%%%%%%%%%%%%%%%%%%%%%
% TimeStamps=TimeStamps-2;
    TimeStamps2=TimeStamps(2:4:length(TimeStamps)); 
    if length(TimeStamps2)<=(trst+trnum-1)
        trnum = length(TimeStamps2)-trst+1;
    end
    trend=trst+trnum-1;
 
    DataTime = (TimeStamps2(2) - TimeStamps2(1));

    cut_spikes = seperate_trials(Spikes,TimeStamps2(trst:trend));    

    BinningTime = [ 0 : BinningInterval : DataTime];
    
        %%%%%% a3 %%%%%%%%%
    x1 = 0:BinningInterval:DataTime-BinningInterval;
    y1=zeros(1,length(x1)); 
    y1(1:2/BinningInterval)=0.18; %(unit:nA)
    y1(4/BinningInterval:6/BinningInterval)=0.18; %(unit:nA)
    y1(8/BinningInterval:10/BinningInterval)=0.18; %(unit:nA)
        
    y1=a_data(3,TimeStamps(2)*20000:TimeStamps(5)*20000);
    x1=1/20000:1/20000:length(y1)/20000;
    %%%% pick diode's timestamps %%%%%% 
%     [~,locs_a2]=findpeaks(diff(a2),'MINPEAKHEIGHT',5*std(diff(a2)));
%     analog_loc = (locs_a2)/1000; 
%     TimeStamps_a2 = analog_loc;
%    
%     [b,a] = butter(2,10/1000,'low'); % set butter filter
%     callumin_filter = filter(b,a,callumin);
%     
%     x1 = 0.004:0.001:DataTime;
%     y1=callumin_filter(TimeStamps_a2(1)*1000:TimeStamps_a2(4)*1000)';    
    
%%%%%%%%%%%%%%%%%%%%%%%%%  Plot Different Trials   %%%%%%%%%%%%%%%%% 
%     if length(TimeStamps2)<=18
%         sweepend=length(TimeStamps2);
%     else
%         sweepend=24;
%     end
    sweepend=trend;
    figure(2);hold on
    set(gcf,'position',[150,30,1024,900])
    h = subplot(sweepend+1,1,1);
    
    for sweepindex=1:sweepend-1
        TimeStampsSweep=TimeStamps2(sweepindex:sweepindex+1); % forcus on which trails 
        cut_spikes = seperate_trials(Spikes,TimeStampsSweep);      
        for i = 1:60
            [n,xout] = hist(cut_spikes{i},BinningTime) ;
            BinningSpike(sweepindex,i,:) = n ;
        end
        subplot(sweepend+1,1,sweepindex);hold on
%         plot(BinningTime,squeeze(sum(BinningSpike(sweepindex,ch,:),2)));
        plot(BinningTime,squeeze(sum(BinningSpike(sweepindex,roi,:),2)));
    end     
        
    subplot(sweepend+1,1,sweepindex+1);
    plot(x1,y1,'r');
    samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
%     ylim([min(y1)-0.01,max(y1)+0.01]);
    set(get(h,'title'),'string',[name,'  ch',num2str(roi)]);
    
    saveas(gcf,[name,'_trials.jpg'],'jpg');
    saveas(gcf,[name,'_trials.fig'],'fig');  

 %%%%%%%%%%%%%%%%% raster plot %%%%%%%%%%%%%%%%%%%
    BinningSpike2 = squeeze(sum(BinningSpike(trst:trend-1,:,:),1));
    figure;hold on
    imagesc(BinningTime,[1:60],BinningSpike2);
    title([name,'(sum over ',sprintf('%.0f',length(TimeStamps2)),' trails) ']);  
    xlabel('t(s)');ylabel('channel');
    colorbar;
    saveas(gcf,[name,'_raster.jpg'],'jpg')
    saveas(gcf,[name,'_raster.fig'],'fig')    
     
%%%%%%%%%%%%%  get peaks  %%%%%%%%%%%%%%%%%%%%%   
%     SB1=sum(BinningSpike(chst1:chend1,:),1)/trnum;  
%     [spks,slocs]=findpeaks(SB1,'minpeakdistance',floor(DataTime/BinningInterval/f),'MINPEAKHEIGHT',p);    
% 
%     if isempty(slocs)
%         slocstime=NaN;
%     end
%     slocstime=(slocs-1)*BinningInterval; 
     
%%%%%%%%%%%%%%%%%%%%%%%%%  Plot histogram   %%%%%%%%%%%%%%%%%        
    figure(3)
    subplot(6,1,2:5); 
    plot(BinningTime,sum(BinningSpike2(roi,:),1)/length(roi)); 
%     title([name,'   Ch',num2str(chst1),' to ',num2str(chend1),sprintf('\n'),' PeakTime(s)=',sprintf('%8.3f',slocstime),sprintf('\n')]);
%     ylabel('firing rate per 5ms');
    subplot(6,1,6);
    plot(x1,y1,'r');
    samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
    xlabel('t(s)');
   
%============================= save figures ===============================    
%     saveas(gcf,[name,'_hist.jpg'],'jpg')
%     saveas(gcf,[name,'_hist.fig'],'fig')
%==========================================================================

%     save(filename,'x','y','chst1','chend1','slocstime', '-append'); % save data into original dat
%     catch
%         [msgstr,msgerr] = lasterr;
%         disp([msgstr,msgerr])
%     end
end


%%  plot single channel
channels=30;
figure('Name',['Ch ',num2str(channels)],'NumberTitle','off')
for ch=channels
    subplot(6,1,2:5); 
    plot(BinningTime,BinningSpike2(ch,:),'LineWidth',1.5); hold on
    title([name,'   Ch ',num2str(channels)]);
    ylabel(['firing rate for bin  ',num2str(BinningInterval),'ms']); 
    xlabel('t(s)');
end
subplot(6,1,6);
plot(x1,y1,'r','LineWidth',1.5);
samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center


% %% plot average firing rate of several channels
% A = [];
% B = [];
% A = []; %selected channels
% B = [];
% data1=[];
% data2=[];
% for i=1:length(A)
%     data1=[data1;BinningSpike2(A(i),:)];
% end
% for i=1:length(B)
%     data2=[data2;BinningSpike2(B(i),:)];
% end
% mean_fr1=sum(data1,1)/length(A);
% mean_fr2=sum(data2,1)/length(B);
% subplot(6,1,2:5); 
% plot(BinningTime,mean_fr1,'Color','r');hold on;
% plot(BinningTime,mean_fr2,'Color','b')
% title([name,'   Ch',num2str(i)]);
% ylabel('firing rate per 5ms');
% subplot(6,1,6);
% plot(x1,y1,'r');
% samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
% xlabel('t(s)');
% 
