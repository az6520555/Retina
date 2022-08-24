% auto correlation of ganglion cell output
clear all
close all

path=['F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200419'];
cd(path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
SamplingRate=20000;
cc=hsv(n_file);

        rr =[9,17,25,33,41,49,...
          2,10,18,26,34,42,50,58,...
          3,11,19,27,35,43,51,59,...
          4,12,20,28,36,44,52,60,...
          5,13,21,29,37,45,53,61,...
          6,14,22,30,38,46,54,62,...
          7,15,23,31,39,47,55,63,...
            16,24,32,40,48,56];
roi = 1:60;
BinningInterval = 0.01;
files=[]; % 

sti_save=[];
fr_save=[];

for z = 1:length(files)
    clearvars -except SamplingRate files BinningInterval roi rr z all_file sti_save fr_save
    file = all_file(files(z)).name 
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    