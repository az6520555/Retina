
clearbin_pos
close all
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210720\merge')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
np=50; %number of repeat of same stimuli period
file_list=[25 24 23 22 21 20 15 14 17];
xparam=[2 1 0.5 0.3 0.1 0.08 0.04 0.025 0]; % background internsity
spike_num=zeros(length(file_list),60);

for z=1:length(file_list)
    file = all_file(file_list(z)).name;
    check_allfilename{z}=file;
    load(file)
    
    for i=1:60
        spike_num(z,i)=length(reconstruct_spikes{i}(reconstruct_spikes{i}>TimeStamps(1) & reconstruct_spikes{i}<TimeStamps(2)));
    end
   
end
figure(1);hold on
figure(2);hold on
for i=1:60
    if max(spike_num(:,i))>300
        figure(1);plot(xparam,spike_num(:,i)/300,'-o')
        figure(2);plot(xparam,spike_num(:,i)/300-spike_num(1,i)/300,'-o')
    end
end
figure(1);
xlabel('background (mW/m^2)')
ylabel('\gamma (Hz)')
figure(2)
xlabel('background (mW/m^2)')
ylabel('\gamma-\gamma_{mean=0.2} (Hz)')
