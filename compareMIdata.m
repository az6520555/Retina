clear all
close all

path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200408\MIdata';
cd(path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
cc=hsv(n_file);

        rr =[9,17,25,33,41,49,...
          2,10,18,26,34,42,50,58,...
          3,11,19,27,35,43,51,59,...
          4,12,20,28,36,44,52,60,...
          5,13,21,29,37,45,53,61,...
          6,14,22,30,38,46,54,62,...
          7,15,23,31,39,47,55,63,... 
            16,24,32,40,48,56];
roi = [1:60];

filelist1=[9 1 7 5 3]; 
for z = 1:length(filelist1)
    file = all_file(filelist1(z)).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    figure(1);subplot(2,1,1)
    hold on;box on
    channel=1;
    title(['channel ',num2str(channel)])
    plot(TimeShift,MI{channel},'LineWidth',2,'LineStyle','-')
    xlim([-1000 1000])
    xlabel('\deltat (ms)')
    ylabel('MI (bits/s)')
    ax = gca;
    ax.XGrid = 'on';
    xline(0)
end
filelist2=[9 1 7 5 3]; 
for z = 1:length(filelist2) 
    file = all_file(filelist2(z)).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    figure(1);subplot(2,1,2)
    hold on;box on
    channel=1;
    title(['channel ',num2str(channel)])
    plot(TimeShift,MI{channel},'LineWidth',2,'LineStyle','-')
    xlim([-1000 1000])
    xlabel('\deltat (ms)')
    ylabel('MI (bits/s)')
    ax = gca;
    ax.XGrid = 'on';
    xline(0)
end


figurepath='F:\§Úªº¶³ºÝµwºÐ\Retina exp\2020 TPS\figure\';
figurename='';
% saveas(gcf,[figurepath,figurename,'.png'])
% saveas(gcf,[figurepath,figurename,'.fig'])