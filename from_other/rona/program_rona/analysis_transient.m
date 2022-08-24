% Analysis transient experiment, by Rona
clear all
close all
clc

layout = [21,19,16,15,12,10,24,22,...
          20,17,14,11, 9, 7,26,25,...
          23,18,13, 8, 6, 5,29,30,...
          28,27, 4, 3, 1, 2,32,31,...
          33,34,57,58,60,59,35,36,...
          38,43,48,53,55,56,37,39,...
          41,44,47,50,52,54,40,42,...
          45,46,49,51];

cd('\\192.168.0.100\Experiment\Public\20180426') ; % the folder of the files
all_file = dir('*.mcd') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file) ; 

for m = 7
    clearvars -except all_file n_file m layout
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    
file = strcat([filename]);
AllDataInfo =datastrm(file); % open data
n_channel = getfield(AllDataInfo,'NChannels2'); % # of channels
n_channel = n_channel(1);% # of channels 

t1=getfield(AllDataInfo,'sweepStartTime');
t2=getfield(AllDataInfo,'sweepStopTime') ;
startend = [t1 t2]; % startend time of all
SamplingRate = getfield(AllDataInfo,'MillisamplesPerSecond2')/1000; % rescale sample rate to Hz
SamplingRate = SamplingRate(2);
TriggerAllData = nextdata(AllDataInfo,'startend',startend,'originorder','on' ,'streamname','Analog Raw Data');

StartStopVector = [t1 t2];
AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data'); % call data, the 60 channels
rawdata = AllData.data ;
[vdata, tvals] = ad2muvolt(AllDataInfo, rawdata, 'Electrode Raw Data'); % rescale vdata(uV) , tvals(us)
tvals = tvals/1E6; % rescale unit to s
% vdata2D = vdata;
% tvals2D = tvals;%reshape(tvals,length(rawdata)/n_channel,n_channel) ;
% Time = tvals2D;

SamplingRate = 20000;
a_data2 = TriggerAllData.data(2,:); 
isi = TriggerAllData.data(3,:); 
BinningInterval=0.010;
SpikeTime = cell(60,1);

[~,locs1]=findpeaks(diff(a_data2),'MINPEAKHEIGHT',5*std(diff(a_data2)));
TimeStamps_h = (locs1)/SamplingRate;
[~,locs2]=findpeaks(-diff(a_data2),'MINPEAKHEIGHT',5*std(diff(a_data2)));
TimeStamps_n = (locs2)/SamplingRate;

H=[];N=[];vdata_H=[];vdata_N=[];
for i = 1 :length(TimeStamps_h)
    H = [H isi(locs1(i):locs2(i))];
    vdata_H = [vdata_H vdata(:,locs1(i):locs2(i))];
    if i == length(TimeStamps_h)
    break;end
    N = [N isi(locs2(i):locs1(i+1))];
    vdata_N = [vdata_N vdata(:,locs2(i):locs1(i+1))];
end

for channel_index = 1:60;
    data = vdata_H(layout(channel_index),:)'; 
    [b,a] = butter(2,200/20000,'high'); % set butter filter
    FilterData = filter(b,a,data);
    [pks,locs] = findpeaks(-FilterData,'MINPEAKHEIGHT', std(FilterData)*4);
    SpikeTime_H{channel_index} = [SpikeTime{channel_index},locs'/SamplingRate] ;
end

for channel_index = 1:60;
    data = vdata_H(layout(channel_index),:)'; 
    [b,a] = butter(2,200/20000,'high'); % set butter filter
    FilterData = filter(b,a,data);
    [pks,locs] = findpeaks(-FilterData,'MINPEAKHEIGHT', std(FilterData)*4);
    SpikeTime_N{channel_index} = [SpikeTime{channel_index},locs'/SamplingRate] ;
end

save(['E:\variables_',name,'.mat'],'vdata','TriggerAllData','H','N','vdata_N','vdata_H','SpikeTime_H','SpikeTime_N','-v7.3')
end