% on off threshold finding
clear
close all
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210720\SplitData')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
np=50; %number of repeat of same stimuli period
loop_ind=[6 8 7 5 9];
contrast=[0.05 0.1 0.15 0.2 0.3];
spike_num=zeros(length(contrast),60);

for z=1:length(loop_ind)
    file = all_file(loop_ind(z)).name;
    check_allfilename{z}=file;
    load(file)
    
    for i=1:60
        spike_num(z,i)=length(Spikes{i}(Spikes{i}>TimeStamps(1) & Spikes{i}<TimeStamps(2)));
    end
   
end
figure;hold on
for i=1:60
    if max(spike_num(:,i))>100
        plot(contrast,spike_num(:,i)/300,'-o')
    end
end
xlabel('contrast')
ylabel('firing rate (Hz)')