clear
filename_mat_path='\\192.168.0.102\Public\Retina\Chou\stimulus_saving\01-Sep-2022';
cd(filename_mat_path)
load(fullfile('stimuli_seq_01-Sep-2022_17-0.mat'))
% create txt file
fid=fopen('20220901_list.txt','wt');

for i=1:size(parameters,1)
    [path,filenames,exts]=fileparts(parameters(i,2));
    fprintf(fid,'%s \n',[filenames,exts]);
end

fclose(fid);
