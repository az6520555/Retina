clear;close all
path='F:\§Úªº¶³ºÝµwºÐ\Master Thesis\Figures\Result\data';
cd(path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);

files=[3 1 4 2];
c=['k','r'];
color_num=[1,2,1,2];
ax(1)=subplot(7,4,[1,2]);
ax(2)=subplot(7,4,[3,4]);
ax(3)=subplot(7,4,[5,6,9,10]);
ax(4)=subplot(7,4,[7,8,11,12]);
ax(5)=subplot(7,4,[5,6,9,10]+8);
ax(6)=subplot(7,4,[7,8,11,12]+8);
ax(7)=subplot(7,4,[5,6,9,10]+16);
ax(8)=subplot(7,4,[7,8,11,12]+16);

r_filenum=[1,1,2,2];
peaktime=zeros(3,4);

for z = 1:length(files)
    file = all_file(files(z)).name
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]); 
    name(name=='_')='-';
    
    subplot(ax(r_filenum(z)));hold on;box on
    plot(tSingleSti,SingleSti,'linewidth',1,'color',c(color_num(z)))
    if r_filenum(z)==1; 
        xlim([0,3]);
        xline(1.5,'--','color',[0.5 0.5 0.5]);
        title('Gaussian pulse (\sigma_G=0.3 s)','fontsize',12);
        ylabel({'Intensity','(mW/m^2)'},'fontsize',12);
    else
        xlim([1.5,4.5]);
        xline(3,'--','color',[0.5 0.5 0.5]);
        title('Gaussian pulse (\sigma_G=0.6 s)','fontsize',12);
    end
    ylim([5,15])
    
    set(gca,'XTick',[])
    for nn=1:3
        subplot(ax(r_filenum(z)+2*nn));hold on;box on
        plot(tBinning,frN{nn},'linewidth',1,'color',c(color_num(z)))
        if r_filenum(z)==1; 
            xlim([0,3]);
            xline(1.5,'--','color',[0.5 0.5 0.5])
            ylabel({['Cell ',num2str(nn)],'Firing Rate (Hz)'},'fontsize',12)
        else
            xlim([1.5,4.5]);
            xline(3,'--','color',[0.5 0.5 0.5])   
        end
        if nn ~= 3;set(gca,'XTick',[]);
        else; xlabel('t (s)','fontsize',12)
        end
    end
    for j=1:3
        [M,imax(j)]=max(frN{j});
        peaktime(j,z)=tBinning(imax(j));
    end
    
end