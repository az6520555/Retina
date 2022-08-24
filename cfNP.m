close all;
clear all;

channel_number = [9 33]; %2 element only for now. 1st for N, 2nd for P.
tolerance = 0.05; %s

load() % load spike time file
load() % load original mutual information
analyze_spikes = sorted_spikes;

%% plot original MI
figure(1);
for k = channel_number
    plot(time,Mutual_infos{k},'LineWidth',1.5,'LineStyle','-');hold on;
    %plot(time,smooth(Mutual_shuffle_infos{j}),'LineWidth',1.5,'LineStyle','-');hold on;
    xlim([ -1500 1500])
    ylim([0 inf+0.5])
    xlabel('\deltat (ms)');ylabel('MI (bits/second)');
    set(gca,'fontsize',12); hold on
    title(strrep(name,'_',' '));
end 



%% remove sharing spike and plot the cancelled, sharing, and original spike
for j = 1:length(channel_number)
    sub_Spikes{j} = analyze_spikes{channel_number(j)};
    if isempty(sub_Spikes{j})==1
        sub_Spikes{j}=0;
    end
    sub_Spikes{j} = sub_Spikes{j}';
    
end
null_index_1 = [];
null_index_2 = [];
%There may be a better way to search
for m = 1:length(sub_Spikes{1})
    for n = 1:length(sub_Spikes{2})
        if abs(sub_Spikes{1}(m)-sub_Spikes{2}(n)) < tolerance
            null_index_1 = [null_index_1 m];
            null_index_2 = [null_index_2 n];
        end
    end
end
null_index_1 = unique(null_index_1);
null_index_2 = unique(null_index_2);
sub_Spikes{3} = sub_Spikes{1}(null_index_1);
sub_Spikes{4} = sub_Spikes{2}(null_index_2);
sub_Spikes{5} = sub_Spikes{1};
sub_Spikes{6} = sub_Spikes{2};
sub_Spikes{1}(null_index_1) = [];
sub_Spikes{2}(null_index_2) = [];
figure;
LineFormat.Color = [0.3 0.3 0.3];
plotSpikeRaster(sub_Spikes,'PlotType','vertline','LineFormat',LineFormat);
hold on;
oN = length(analyze_spikes{channel_number(1)});
oP = length(analyze_spikes{channel_number(2)});
analyze_spikes{channel_number(1)}(null_index_1) = [];
analyze_spikes{channel_number(2)}(null_index_2) = [];
cN = length(analyze_spikes{channel_number(1)});
cP = length(analyze_spikes{channel_number(2)});
cN/oN
cP/oP

%% calculate and plot  MI of cancelled spike

TheStimuli=bin_pos;
bin=BinningInterval*10^3; %ms
BinningTime =diode_BT;

StimuSN=30; %number of stimulus states
nX=sort(TheStimuli);
abin=length(nX)/StimuSN;
intervals=[nX(1:abin:end) inf]; % inf: the last term: for all rested values
temp=0; isi2=[];
for jj=1:length(TheStimuli)
    temp=temp+1;
    isi2(temp) = find(TheStimuli(jj)<=intervals,1);
end

    %% BinningSpike
BinningSpike = zeros(60,length(BinningTime));
for  j = channel_number  % i is the channel number
    [n,~] = hist(analyze_spikes{j},BinningTime) ;  %yk_spikes is the spike train made from"Merge_rePos_spikes"
    BinningSpike(j,:) = n ;
end

    %% Predictive information
backward=ceil(15000/bin);
forward=ceil(15000/bin);
time=[-backward*bin:bin:forward*bin];
figure(1);
for j = channel_number
    
    Neurons = BinningSpike(j,:);  %for single channel
    information = MIfunc(Neurons,isi2,BinningInterval,backward,forward);
    plot(time,information,'LineWidth',1.5,'LineStyle','-');hold on;
    %plot(time,smooth(Mutual_shuffle_infos{j}),'LineWidth',1.5,'LineStyle','-');hold on;
    xlim([ -1500 1500])
    ylim([0 inf+0.5])
    set(gca,'fontsize',12); hold on
    %         legend('-DynamicLegend');
    %         legend('show')
    %lgd = legend(['Original N'],['Original P'], ['Cancelled N'],['Cancelled P']);
end

analyze_spikes{channel_number(1)} = sorted_spikes{channel_number(1)}(null_index_1);
analyze_spikes{channel_number(2)} = sorted_spikes{channel_number(2)}(null_index_1);


TheStimuli=bin_pos;
bin=BinningInterval*10^3; %ms
BinningTime =diode_BT;

StimuSN=30; %number of stimulus states
nX=sort(TheStimuli);
abin=length(nX)/StimuSN;
intervals=[nX(1:abin:end) inf]; % inf: the last term: for all rested values
temp=0; isi2=[];
for jj=1:length(TheStimuli)
    temp=temp+1;
    isi2(temp) = find(TheStimuli(jj)<=intervals,1);
end

    %% BinningSpike
BinningSpike = zeros(60,length(BinningTime));
for  j = channel_number  % i is the channel number
    [n,~] = hist(analyze_spikes{j},BinningTime) ;  %yk_spikes is the spike train made from"Merge_rePos_spikes"
    BinningSpike(j,:) = n ;
end

    %% Predictive information
backward=ceil(15000/bin);
forward=ceil(15000/bin);
time=[-backward*bin:bin:forward*bin];
figure(1);
for j = channel_number
    
    Neurons = BinningSpike(j,:);  %for single channel
    information = MIfunc(Neurons,isi2,BinningInterval,backward,forward);
    plot(time,information,'LineWidth',1.5,'LineStyle','-');hold on;
    %plot(time,smooth(Mutual_shuffle_infos{j}),'LineWidth',1.5,'LineStyle','-');hold on;
    xlim([ -1500 1500])
    ylim([0 inf+0.5])
    set(gca,'fontsize',12); hold on
    lgd = legend(['Original N'],['Original P'], ['Cancelled N'],['Cancelled P'], ['sharing N'],['sharing P']);
end
    
    
    
    
    
    
    
    
    
