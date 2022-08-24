
clear
newname_path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\temporal_pattern\20200508\08-May-2020';
cd(newname_path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
for i =1:n_file;date{i,1}=all_file(i).date;end
[datesort ind]=sort(date)

oldname_path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\temporal_pattern\20200508\SplitData'
cd(oldname_path)
for i=1:n_file
    newname=all_file(ind(i)).name
    movefile([num2str(i),'.mat'],newname(7:end))
end