clear all
% close all

%%%%%%%%%%%%%%%%%%%%%% user's setting %%%%%%%%%%%%%%%%%%%%%%%%%%%

BinningInterval = 5/1000 ; % (s) binning duration = 5ms
n_channel = 60;

%%%%%%%%%%%%%%%%%%%%%% call all data in %%%%%%%%%%%%%%%%%%%%%%%%%% 
cd('E:\Chou\20210504') ; % the folder of the files
all_file = subdir('*.mcd') ; 
n_file = length(all_file) ; 
files=1:n_file;

for m = [1]
    clearvars -except all_file n_file m BinningInterval SamplingRate 
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file)
    directory = [pathstr,'\']
    filename = [name,ext]
    
%     cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\analysis code\rona')
%     [Spikes,TimeStamps,a_data,Infos] = analyze_MEA_data([directory,filename],1,'','Chou','all',100000); %filename,save_data,comment,experimenter,analog_type,r_interval
    [reconstruct_spikes,TimeStamps,a_data,Infos] = analyze_MEA_data_revised([directory,filename],1,'','Chou','all',100000);
    
% a_data = ad2muvolt(AllDataInfo,a_data,'Electrode Raw Data');%a_data2=(a_data-32768)*0.1042;
    
% %%%%%%%%%  belows are written in  "seperate_trials" "analyze_MEA_data"  %%%%%%%%%
%     file = strcat([directory,filename]);
%     AllDataInfo =datastrm(file); % open data
%     n_channel = getfield(AllDataInfo,'NChannels2'); % # of channels
%     n_channel = n_channel(1);% # of channels 
% 
%     t1=getfield(AllDataInfo,'sweepStartTime');
%     t2=getfield(AllDataInfo,'sweepStopTime') ;
%     startend = [t1 t2]; % startend time of all
%     SamplingRate = getfield(AllDataInfo,'MillisamplesPerSecond2')/1000; % rescale sample rate to Hz
%     SamplingRate = SamplingRate(2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% call trigger time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     TriggerAllData = nextdata(AllDataInfo,'startend',startend,'originorder','on' ,'streamname','Analog Raw Data');
%     TriggerData = TriggerAllData.data(1,:); % TriggerAllData.data(1,:) for Kevin   TriggerAllData.data(1,:) for Micheal
%     [pks,locs] = findpeaks(diff(TriggerData),'MINPEAKHEIGHT',100);  % get the trigger time (locs)
%     clearvars pks locs 
% 
%     read_interval = 1000000; 
%     SpikeTime = cell(60,1);
% 
% for i = t1 : read_interval : t2-1 %% get all data into one vector
%     if i + read_interval-1 < t2-1
%         StartStopVector = [i , i+read_interval-1];
%     else
%         StartStopVector = [i , t2-1];
%     end
%     
%     %%%%%  for data before 2015/9/14 %%%%%%%%%%%
% %     AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data 1'); % call data, the 60 channels
% %     rawdata = AllData.data ;
% %     [vdata, tvals] = ad2muvolt(AllDataInfo, rawdata', 'Electrode Raw Data 1'); % rescale vdata(uV) , tvals(us)
% %     tvals = tvals/1E6; % rescale unit to s
% %     vdata2D = reshape(vdata,n_channel,length(rawdata)/n_channel) ; % separate 60 channels
% %     tvals2D = reshape(tvals,length(rawdata)/n_channel,n_channel) ;
% %     Time = tvals2D(:,1)';
% %     
% %%%%%%  for data after 2015/9/14 %%%%%%%%%%%
%     AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data'); % call data, the 60 channels
%     rawdata = AllData.data ;
%     [vdata, tvals] = ad2muvolt(AllDataInfo, rawdata, 'Electrode Raw Data'); % rescale vdata(uV) , tvals(us)
%     tvals = tvals/1E6; % rescale unit to s
%     vdata2D = vdata;
%     tvals2D = tvals;%reshape(tvals,length(rawdata)/n_channel,n_channel) ;
%     Time = tvals2D;
%    
% 
% %%%%%%%%%%%%%%%%%%%%%%%% process channel by channel %%%%%%%%%%%%%%%%%%%%%
%    
%     BinningTime = [0:BinningInterval:Time(end)] ;
%     for channel_index = 1:n_channel;
%             data = vdata2D(channel_index,:)';
%             [b,a] = butter(2,200/10000,'high'); % set butter filter
%             FilterData = filter(b,a,data);
%             %     figure(2),subplot(211),plot(Time,data),subplot(212),plot(Time,FilterData)
%             std_data = std(FilterData);
%             [pks,locs] = findpeaks(-FilterData,'MINPEAKHEIGHT', std_data*5);
%             locs = locs+i/1000*SamplingRate;
%             %     plot(-FilterData), hold on, plot(locs,std_data*5,'o'),hold off
%             SpikeTime{channel_index} = [SpikeTime{channel_index},locs'/SamplingRate] ;
%     end
% 
% end
%%%%%%%%%%%%%%   seperate trails    %%%%%%%%%%%%%%%%%%%%%%
% for i = 1:n_channel
%     temp = Spikes{i};
%     temp(temp<TimeStamps(1)) = [];
%     for j = 1:length(temp)
%         loc = find(temp(j)<TimeStamps,1)-1;
%         if isempty(loc) || loc == 0
%             continue
%         end
%         temp(j) = temp(j) - TimeStamps(loc);
%     end
%     temp(temp>TimeStamps(end))=[];
%     temp(temp==0)=[];
%     SpikeTime2{i} = temp;
% end
% % 
% 
% spike_rearrange = cell(size(SpikeTime2));
% for i = 1:60
%     spike_rearrange{i} = SpikeTime2{MEA_layout_Tina03(i)};
% end
% DataTime = (TimeStamps(2)-TimeStamps(1));
% 
% BinningTime = [0:BinningInterval:DataTime] ;
% BinningSpike = zeros(60,length(BinningTime)) ;
% 
% for i = 1:60
%     [n,xout] = hist(spike_rearrange{i},BinningTime) ;
%     BinningSpike(i,:) = n ;
% end



% x = 0:BinningInterval:DataTime;
% y = 100*((sin(2*pi*x/DataTime)+1)/2+2); %% M=2.2V;m=1.2V

% [p,t]=max(sum(BinningSpike));
% wf=(y(t+1)-y(t-1))/y(t);  %%weber's fraction


% 
% figure
% plot(BinningTime,cut_spikes/length(TimeStamps));
% % xlim([0,2])
% % title(['Continuous Sinewave Stimulus(3mins)','   ','T=', num2str(DataTime),'s','   ','WF=', num2str(wf)]);
% title([char(filename)]);
% xlabel('t(s)');ylabel('spike# with 60 channels in 2ms');
% saveas(gcf,[filename(1:end-4),'_hist.jpg'],'jpg')
% saveas(gcf,[filename(1:end-4),'_hist.fig'],'fig')
% 
% 
% % % 
% 
% figure(1)
% hold on
% plot(BinningTime/BinningTime(end),sum(BinningSpike),'color',[rand(1),rand(1),rand(1)] ,'DisplayName', [filename(1:end-4)]); % normalize
% % legend([filename(1:end-4)]);
end
% x = 0:BinningInterval:1;
% y = 10*((sin(2*pi*x)+1)/2+1.2);
% % plot(x,y)
% legend(gca,'show')
% hold off

    
   





