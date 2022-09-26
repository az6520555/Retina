%  LN model calculate nonlinearity
clear 
close all
datapath='\\192.168.0.102\Public\Retina\Chou\Exp\20220916\SplitData';
cd(datapath);
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
SamplingRate = 20000;
rr =    [9,17,25,33,41,49,...
      2,10,18,26,34,42,50,58,...
      3,11,19,27,35,43,51,59,...
      4,12,20,28,36,44,52,60,...
      5,13,21,29,37,45,53,61,...
      6,14,22,30,38,46,54,62,...
      7,15,23,31,39,47,55,63,...
        16,24,32,40,48,56];
roi = [1:60];

file_numbers=[16 18 17 19 20];
colors_default={[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560], ...
    [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};

for z=1:length(file_numbers)
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr...
        information rr STA roi fignames date datapath Date colors_default file_numbers ls_set i_set
    
    filename=all_file(file_numbers(z)).name;
    load(filename);
    bin = 10;  BinningInterval = bin*10^-3;  %ms
    %% stimulus
    [b,a] = butter(2,50/20000,'low'); % set butter filter
    a_data2 = filter(b,a,a_data(1,:));
    TriggerData = a_data2(TimeStamps(1)*SamplingRate:TimeStamps(length(TimeStamps))*SamplingRate);% figure;plot(isi);
    inten = downsample(TriggerData,SamplingRate*BinningInterval);









end
