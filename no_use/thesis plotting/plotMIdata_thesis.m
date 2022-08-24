clear all
close all

path='F:\我的雲端硬碟\Retina exp\exp data\Sorted_final_data\20200408\MIandSTA';
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
filelistMI=[23 18 13 8];  
filelistSTA=[];
MI_savedata=[];

for z = 1:length(filelistMI)
    file = all_file(filelistMI(z)).name;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    TS=TimeShift;
    
    %% plot multichannel data
%     figure(1);
%     for i=1:60
%         subplot(8,8,rr(i));hold on
%         plot(TimeShift,MI{i},'LineWidth',1,'LineStyle','-');
%         xlim([-1000 1000])
%     end
%     
%     file_onoff='F:\我的雲端硬碟\Retina exp\exp data\整理\OnOff_index\20200318-onoff2-sort-unit2onoff_index.mat';
%     onoff_color(file_onoff)
%     file_NP='F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\20200318_OU_tau=op5_2_sort_unit2_MI_also_other_parameters.mat'
%     NPcolor(file_NP)
    
    %% plot singlechannel data
    figure(3);hold on;box on
    channel=[53,36,32];
%     title(['channel ',num2str(channel)])
    for ii=1:3
        ax(ii)=subplot(1,6,[ii*2-1,ii*2]);hold on;box on
        plot(TimeShift,MI{channel(ii)},'LineWidth',1.5,'LineStyle','-')
        xlim([-1000 1000])
        xlabel('\deltat (ms)')
        ax = gca;
        if ii ==1;
            ylabel('MI (bits/s)')
        end
%         if ii==3
%             legend('OU not filtered ','f_c=3.5 Hz','f_c=2 Hz','f_c=1 Hz')
%         end
    %     ax.XGrid = 'on';
        
    end
    
    %% plot figure of every channel and STA
%     file = all_file(filelistSTA(z)).name ;
%     [pathstr, name, ext] = fileparts(file);
%     directory = [pathstr,'\'];
%     filename = [name,ext];
%     load([filename]);
%     mkdir MIfigure
%     for channel=1:60
%         figure(channel)
% %         subplot(2,1,1);hold on;box on
%         plot(TS,MI{channel},'LineWidth',2,'LineStyle','-')
%         xlabel('\deltat (ms)')
%         ylabel('MI (bits/s)')
%         title(['MI  channel ',num2str(channel)])
%         ax = gca;ax.XGrid = 'on';xline(0)
%         subplot(2,1,2);hold on;box on
%         plot(TimeShift{channel},STAAAAA{channel},'LineWidth',2,'LineStyle','-')
%         title(['STA  channel ',num2str(channel)])
%         xlim([-1000 1000])
%         xlabel('time to spike (ms)')
%         ylabel('STA')
end
for i =1:3    
    subplot(subplot(1,6,[i*2-1,i*2]))
    xline(0,'--','color',[0.5,0.5,0.5])
end
legend('OU not filtered ','f_c=3.5 Hz','f_c=2 Hz','f_c=1 Hz')