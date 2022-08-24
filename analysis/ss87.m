clear
close all
path='\\ZebraNas\Public\Retina\troy\SplitData';
cd(path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
for z=[6 7 8 9]
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    plot(a_data(1,1:10:end));hold on
end