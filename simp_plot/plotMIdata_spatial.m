clear all
close all

path='\\192.168.0.102\Public\Retina\Chou\Exp\data_2021\20210506\MI\unsort';
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
  
% NPchannel_date='20210506';
roi = [1:60];
filelistMI=[]; 
filelistSTA=[];
MI_savedata=[];


for z = 1:length(filelistMI)
    file = all_file(filelistMI(z)).name;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    TS=time;
    
    %% plot multichannel data
    f1=figure(1);
    f1.Position=[10 50 1900 950]
    if z==1 
        for i=1:60
            subplot(8,8,rr(i));hold on
            xline(0,':') 
        end
    end
    for i=1:60
        subplot(8,8,rr(i));hold on
        plot(TS,Mutual_infos{i},'LineWidth',1,'LineStyle','-','DisplayName',all_file(filelistMI(z)).name(1:end-4));
        xlim([-1000 1000])
        title(num2str(i))
    end
    legend('Interpreter', 'none', 'Position', [0.89 0.92 0.03 0.02])
    

%     file_onoff='F:\我的雲端硬碟\Retina exp\exp data\整理\OnOff_index\20200318-onoff2-sort-unit2onoff_index.mat';
%     onoff_color(file_onoff)
%     path_NP='F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\NP_spatial';
%     NPcolor(fullfile(path_NP,[NPchannel_date,'_NP.mat']))
    
    if z==length(filelistMI)
        mkdir MI_figure
        cd('MI_figure')
        set(gcf, 'InvertHardCopy', 'off');
        print(gcf,[path,'\MI_figure\all_channel.png'],'-dpng','-r300')
    end
    
    %% plot singlechannel data
    figure(3);hold on;box on
    channel=32;
    title(['channel ',num2str(channel)])
    plot(TS,Mutual_infos{channel},'LineWidth',2,'LineStyle','-')
    xlim([-1000 1000])
%     legend('Interpreter', 'none','FontSize',4)
    xlabel('\deltat (ms)')
    ylabel('MI (bits/s)')
    ax = gca;
    ax.XGrid = 'on';
    %% plot figure of every channel 
%     for channel=1:60
%         figure(channel);hold on;box on
%         plot(TS,Mutual_infos{channel},'LineWidth',2,'LineStyle','-')
%         xlabel('\deltat (ms)')
%         ylabel('MI (bits/s)')
%         title(['MI  channel ',num2str(channel)])
%         ax = gca;ax.XGrid = 'on';xline(0)
%         xlim([-1000 1000])
%     end
%     
%     if z==length(filelistMI)
%         mkdir MI_figure
%         cd('MI_figure')
%         for channel=1:60
%             figure(channel)
%             print(gcf,['channel',num2str(channel),'.png'],'-dpng','-r300')
%         end
%     end
        
    %% save data
%     saving_channel=39;
%     MI_savedata=[MI_savedata MI{saving_channel}']; % MI

end

% cd([path,'\MIfigure']);mkdir(name);cd([path,'\MIfigure\',name])
% for i=1:60;figure(i);
%     legend('no filter','fc=10','fc=7','fc=4','fc=2');saveas(gcf,[path,'\MIfigure\',name,'\',name,'_channel',num2str(i),'.png']);
% end


%% save xls file
% MI_savedata=[TS' MI_savedata];
% path_data='F:\我的雲端硬碟\Retina exp\Retina as NGD paper submit\data';
% data_name='0419 data.xls'
% writematrix(MI_savedata,fullfile(path_data,data_name))