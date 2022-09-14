clear
filename_mat_path='\\192.168.0.102\Public\Retina\Chou\Exp\20220901\01-Sep-2022';
cd(filename_mat_path)
all_file=dir('*.mat');
n_file=length(all_file);


% create txt file
fid=fopen('20220901_list.txt','wt');

for i=1:size(parameters,1)
    [path,filenames,exts]=fileparts(parameters(i,2));
    fprintf(fid,'%s\n',[filenames]);
end

fclose(fid);
