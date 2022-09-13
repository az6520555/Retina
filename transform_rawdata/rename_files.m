
clear
newname_path='\\192.168.0.102\Public\Retina\Chou\Exp\20220907\07-Sep-2022'; % change every time
cd(newname_path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
for i =1:n_file;date{i,1}=all_file(i).date;end
[datesort ind]=sort(date)

oldname_path='\\192.168.0.102\Public\Retina\Chou\Exp\20220907\SplitData'  % change every time
cd(oldname_path)
for i=1:n_file
    movefile([num2str(i),'.mat'],all_file(ind(i)).name)
end