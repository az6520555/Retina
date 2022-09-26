%  LN model calculate nonlinearity
clear 
close all
datapath='\\192.168.0.102\Public\Retina\Chou\Exp\20220916\SplitData';
cd(datapath);
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
SamplingRate = 20000;