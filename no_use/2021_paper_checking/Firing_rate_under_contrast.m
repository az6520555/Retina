
clear
close all
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210720\SplitData')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
np=50; %number of repeat of same stimuli period
loop_ind=[6 8 7 5 9]; % select files 
xparam=[0.05 0.1 0.15 0.2 0.3]; % specify the contrasts
spike_num=zeros(length(xparam),60);

for z=1:length(loop_ind)
    file = all_file(loop_ind(z)).name;
    check_allfilename{z}=file;
    load(file)
    
    for i=1:60
        spike_num(z,i)=length(Spikes{i}(Spikes{i}>TimeStamps(1) & Spikes{i}<TimeStamps(2)));
    end
   
end
figure(1);hold on
figure(2);hold on
for i=1:60
    if max(spike_num(:,i))>100
        figure(1);plot(xparam,spike_num(:,i)/300,'-o')
        figure(2);plot(xparam,spike_num(:,i)/300-spike_num(1,i)/300,'-o')
    end
end
figure(1);
xlabel('contrast')
ylabel('\gamma (Hz)')
figure(2)
xlabel('contrast')
ylabel('\gamma-\gamma_{C=0.05} (Hz)')