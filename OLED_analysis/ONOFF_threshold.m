% on off threshold finding
clear
close all
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210504\data')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
np=50; %number of repeat of same stimuli period
loop_ind=[27 28 29 30 31 35 32 33 34 38 36 37 39];
rate=20000;
bin=0.01;
fig=zeros(1,60);
check_allfilename=cell(length(loop_ind),1);
pulse_inten=[0.05 0.1 0.3 0.5 0.7 1 1.3 1.6 1.9 2 2.2 2.5 3];

for z=1:length(loop_ind)
    file = all_file(loop_ind(z)).name;
    check_allfilename{z}=file;
    load(file)
    % ========= firing rate ============
    BinningTime=0:bin:size(a_data,2)/rate;
    [n xout]=hist(Spikes{11},BinningTime);
    
    % ========= find the number of spikes in a pulse ==========
    time=1/rate:1/rate:size(a_data,2)/rate;
    sti=(a_data(3,:)-min(a_data(3,:)));
    sti=sti/max(sti);
    rectified_sti=-abs(sti-0.5);
    [~,half_locs]=findpeaks(rectified_sti,'MINPEAKDISTANCE',100,'MINPEAKHEIGHT',-0.2);
    onlocs=half_locs(1:2:end-1);
    offlocs=half_locs(2:2:end);
    ON_spike_number=zeros(60,length(onlocs));
    OFF_spike_number=zeros(60,length(onlocs));
    for channel=1:60
        for num_pulse=1:length(onlocs)
            ON_spike_time_inpulse=find(Spikes{channel}>time(onlocs(num_pulse)) & Spikes{channel}<time(offlocs(num_pulse)));
            ON_spike_number(channel,num_pulse)=length(ON_spike_time_inpulse);
        end
        for num_pulse=1:length(offlocs)-1
            OFF_spike_time_inpulse=find(Spikes{channel}>time(offlocs(num_pulse)) & Spikes{channel}<time(onlocs(num_pulse+1)));
            OFF_spike_number(channel,num_pulse)=length(OFF_spike_time_inpulse);
        end
        ON_firing(channel,z)=sum(ON_spike_number(channel,5:end))/length(loop_ind);
        OFF_firing(channel,z)=sum(OFF_spike_number(channel,5:end))/length(loop_ind);
    end
%     figure;plot(ON_spike_number);hold on;plot(OFF_spike_number)
%     figure;yyaxis left;plot(xout,n);yyaxis right;plot(time,sti)
end
for ch=1:60
    figure(1);plot(pulse_inten,ON_firing(ch,:),'-o');hold on
    figure(2);plot(pulse_inten,OFF_firing(ch,:),'-o');hold on
end