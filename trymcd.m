clear
close all
rate=20000;
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210420\data\raw')
filename='0224_Period_doubling_10Hz.mcd'; %20200508_merge
file = strcat(filename);
AllDataInfo =datastrm(file);
SelectedChannel=1:60;

layout = [21,19,16,15,12,10,24,22,...
          20,17,14,11, 9, 7,26,25,...
          23,18,13, 8, 6, 5,29,30,...
          28,27, 4, 3, 1, 2,32,31,...
          33,34,57,58,60,59,35,36,...
          38,43,48,53,55,56,37,39,...
          41,44,47,50,52,54,40,42,...
          45,46,49,51];

streamname=getfield(AllDataInfo,'StreamNames')
t1=getfield(AllDataInfo,'sweepStartTime');
t2=getfield(AllDataInfo,'sweepStopTime');
StartStopVector=[t1 t2];

%% ======================== electric raw data ======================================
% get stimulation form
AllData2 = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Analog Raw Data');
adata = AllData2.data;
adatasti=adata(3,:);
datatime=0:1/rate:(length(adatasti)-1)/rate;
adatasti=adatasti(1:10:end);
stitime=datatime(1:10:end);
% get response raw data
AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data'); % call data, the 60 channels
rawdata = AllData.data;
[vdata, ~] = ad2muvolt(AllDataInfo, rawdata, 'Electrode Raw Data');
for num=1:60
    data=vdata(layout(SelectedChannel(num)),:);
    [b,a] = butter(2,[1 30]/(rate/2),'bandpass'); % set butter filter
    FilterData = filter(b,a,data);

    %% separate repeat trials
    % beginind=10*rate;
    % adata_kai=-adata(1,beginind:end);
    % adata_kai=adata_kai-min(adata_kai);
    % adatap=findpeaks(adata_kai,'minpeakheight',1/2*max(adata_kai),'minpeakdistance',1000);
    % figure;
    % plot(adata_kai);hold on
    % plot(adatap,zeros(length(adatap)),'o')



    %% Fourier
    % ========= select section of ERG ==============
    SelectedSection=[20 60]; %second
    dataind=rate*SelectedSection(1):rate*SelectedSection(2);

    % X=FilterData(dataind);
    X=data(dataind);

    Fs=rate;
    T=1/Fs;
    L=length(X);
    t=(0:L-1)*T;

    Y=fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    f = Fs*(0:(L/2))/L;
    
    % ======== plot data ============
    figure(SelectedChannel(num));
    subplot(2,1,1)
    title(['filtered data, Channel ',num2str(SelectedChannel(num))])
    yyaxis left
    plot(stitime,adatasti)
    yyaxis right
    plot(datatime,FilterData)
    % plot(datatime,data)
    xlim([20,25])
    
    subplot(2,1,2)
    title(['FFT, Channel ',num2str(SelectedChannel(num))])
    plot(f,P1)
    xlim([0 20])
end