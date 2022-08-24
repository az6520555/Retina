clear all;
code_folder = pwd;
exp_folder_cell = {'F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210413\'};
type_folder_cell = {'pos', 'v', 'pos&v'};%'abs', 'pos', 'v', 'pos&v'.
for ttt = 1
for eee = 1
exp_folder = exp_folder_cell{eee};
cd(exp_folder);
mkdir MI
cd MI
inputtype = type_folder_cell{ttt}; 
type = inputtype;
sorted = 0;
unit = 1; %unit = 0 means using 'unit_a' which is writen down while picking waveform in Analyzed_data.

if sorted
    mkdir sort
    %cd(exp_folder)
    cd ([exp_folder,'\MI\sort'])
    all_file = subdir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
    n_file = length(all_file) ;
else
    mkdir unsort
    cd ([exp_folder,'\merge'])
    all_file = subdir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
    n_file = length(all_file) ;
end
cd(code_folder);


