clear
filename_mat_path='\\192.168.0.102\Public\Retina\Chou\stimulus_saving\01-Sep-2022';
load(fullfile(filename_mat_path,'stimuli_seq_01-Sep-2022_17-0.mat'))

oldname_path='\\192.168.0.102\Public\Retina\Chou\Exp\20220823\SplitData';  % change every time
cd(oldname_path)
all_file=dir('*.mat');
n_file=length(all_file);
if n_file==size(parameters,1)
    for i=1:n_file
        movefile([num2str(i),'.mat'],all_file(ind(i)).name)
    end
else
    'error with number of files';
end