
clear
newname_path='\\192.168.0.100\Public\chou\stimulus_saving\20-Jul-2021'; % change every time
cd(newname_path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
for i =1:n_file;date{i,1}=all_file(i).date;end
[datesort ind]=sort(date)

oldname_path='E:\Chou\20210720\whole_field\SplitData'  % change every time
cd(oldname_path)
for i=1:n_file
    movefile([num2str(i),'.mat'],all_file(ind(i)).name)
end