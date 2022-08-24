%plot LPOU different contrast
clear all
close all

path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210720\SplitData\MIandSTA';
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
filelistMI=[7 11 9 5 13];  
filelistSTA=[16 18 17 19 20];
MI_savedata=[];
colorset={[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250]	,[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880]};

for z = 1:length(filelistMI)
    file = all_file(filelistMI(z)).name;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    TS_MI=TimeShift;
    file = all_file(filelistSTA(z)).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    TS_STA=TimeShift{1};
%     

    
    %% plot OFF
    figure(3);hold on;box on
    channel=[9,17,51];
%     title(['channel ',num2str(channel)])
    for ii=1:3
        subplot(3,2,ii*2-1);hold on;box on
        plot(TS_MI,MI{channel(ii)},'LineWidth',1.5,'LineStyle','-','color',colorset{z})
        xlim([-1000 1000])
        if ii==3;xlabel('\delta t (ms)');end
%         if ii==1;title('OFF cell');end
        ylabel('MI (bits/s)')
        xline(0,'--','color',[0.8 0.8 0.8],'linewidth',1)
        
        subplot(3,2,ii*2);hold on;box on
        plot(TS_STA,STAAAAA{channel(ii)},'LineWidth',1.5,'LineStyle','-','color',colorset{z})
        xlim([-600 0])
        if ii==3;xlabel('time to spike (ms)');end
        ylabel('normalized STA')
        yline(0,'--','color',[0.8 0.8 0.8],'linewidth',1)
    end
    sgtitle('OFF cell','fontweight','bold') 
    set(gcf, 'Position',  [100, 100, 600, 800])
    
%% plot ON
    figure(4);hold on;box on
    channel=[30, 52, 56];
%     title(['channel ',num2str(channel)])
    for ii=1:3
        subplot(3,2,ii*2-1);hold on;box on
        plot(TS_MI,MI{channel(ii)},'LineWidth',1.5,'LineStyle','-','color',colorset{z})
        xlim([-1000 1000])
        if ii==3;xlabel('\delta t (ms)');end
%         if ii==1;title('OFF cell');end
        ylabel('MI (bits/s)')
        xline(0,'--','color',[0.8 0.8 0.8],'linewidth',1)
        
        subplot(3,2,ii*2);hold on;box on
        plot(TS_STA,STAAAAA{channel(ii)},'LineWidth',1.5,'LineStyle','-','color',colorset{z})
        xlim([-600 0])
        if ii==3;xlabel('time to spike (ms)');end
        ylabel('normalized STA')
        yline(0,'--','color',[0.8 0.8 0.8],'linewidth',1)
        
    end
    sgtitle('ON cell','fontweight','bold') 
    set(gcf, 'Position',  [1000, 100, 600, 800])
end
% for i =1:3    
%     subplot(subplot(1,6,[i*2-1,i*2]))
%     xline(0,'--','color',[0.5,0.5,0.5])
% end
% legend('OU not filtered ','f_c=3.5 Hz','f_c=2 Hz','f_c=1 Hz')