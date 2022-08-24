%compare the spike trial of same stimuli
clear
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\temporal pattern\20190419')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
k=0;
loop_ind=[1 2 3];
for z = loop_ind
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA roi fignames date path_sort k
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';

    for i=1:60
        if isempty(Spikes{i})==1
            Spikes{i}=0;
        end
    end
    for j=1:60
        Spikes{j}=Spikes{j}-TimeStamps(1);
        pre_sti=find(Spikes{j}<0);
        Spikes{j}(pre_sti)=0;
    end
    
    k=k+1;
    figure(87);
    subplot(3,1,k)
    LineFormat.Color = [0.3 0.3 0.3];
    plotSpikeRaster(Spikes,'PlotType','vertline','LineFormat',LineFormat);
    hold on
   
end
 samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center