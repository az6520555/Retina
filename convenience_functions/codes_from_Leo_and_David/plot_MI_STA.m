close all;
clear all;
% set(0,'DefaultFigureVisible','off')
code_folder = pwd;
exp_folder = 'F:\�ڪ����ݵw��\Retina exp\exp data\Sorted_final_data\20200419\PCAdata';
% exp_folder = 'D:\Leo\0409';
save_photo = 1;
cd(exp_folder);
load('OU_tau=100ms_cutoff=4_19-Apr-2020_0_sort_unit1.mat')
% load('C:\calibration\20190817\boundary_set.mat')
% load('RGC.mat')
%load('different_G.mat')
% load('predictive_channel\bright_bar.mat')
mkdir STA\FIG
cd STA\MI
forward = 100;%90 bins before spikes for calculating STA
backward = 90;%90 bins after spikes for calculating STA
roi = ;
%roi = 11;
for z = 1:n_file
    %choose file
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load(filename);
    name=[name];
    z
    name
    direction = Get_direction(name);
    if strcmp(direction,'UL_DR') || strcmp(direction,'UR_DL')
        micro_per_pixel_long = micro_per_pixel;
    else
        micro_per_pixel_long = micro_per_pixel;
    end
    for i = roi
        if isempty(positive_PCAs{i}) || isempty(negative_PCAs{i})
            continue
        end
        
        figure('units','normalized','outerposition',[0 0 1 1])
        ha = tight_subplot(1,3,[.04 .04],[0.09 0.02],[.04 .04]);
        set(ha, 'Visible', 'off');
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        fig =gcf;
        fig.PaperPositionMode = 'auto';
        fig.InvertHardcopy = 'off';
        axes(ha(1));
%         if sum(RGCs{i}.center_RF) >0
%             newXpos = Monitor2DCoor2BarCoor(RGCs{i}.center_RF(1),RGCs{i}.center_RF(2),direction,'OLED');
% %             newXpos = Monitor2DCoor2BarCoor(RGCs{i}.center_RF(1),RGCs{i}.center_RF(2),direction,'LCD');
%             plot(STA_time,(mean(positive_PCAs{i},1)-newXpos)*micro_per_pixel_long,'r-.');hold on%Red is positive_PCA
%             plot(STA_time,(mean(negative_PCAs{i},1)-newXpos)*micro_per_pixel_long,'b-.');%Blue is negative_PCA
%         else
%             plot(STA_time,(mean(positive_PCAs{i},1)-meaCenter_x)*micro_per_pixel_long,'r');hold on%Red is positive_PCA
%             plot(STA_time,(mean(negative_PCAs{i},1)-meaCenter_x)*micro_per_pixel_long,'b');%Blue is negative_PCA
%         end
        plot(STA_time,(mean(positive_PCAs{i},1)),'r');hold on%Red is positive_PCA
        plot(STA_time,(mean(negative_PCAs{i},1)),'b');%Blue is negative_PCA
        xlabel('time before spike(sec)')
        ylabel('STA from two kinds of clusters')
        title('STA of two different clusters using k means')
        xline(0)
        axes(ha(2)); 
%         if sum(RGCs{i}.center_RF) >0
%             plot(STA_time,(dis_STA(i,:)-newXpos)*micro_per_pixel_long,'k-.');
%         else
            plot(STA_time,(dis_STA(i,:)));
%         end
        xlabel('time before spike(sec)')
        ylabel('STA of bar')
        xline(0)
%         if ismember(i,p_channel)
%             title(['p',int2str(i)])
%         elseif ismember(i,np_channel)
%             title(['np',int2str(i)])
%         end
        legend([num2str(corr_time),' sec'])
        axes(ha(3));
        plot(time,smooth(B_pos_Mutual_infos{i}),'b');hold on
        plot(time,smooth(R_pos_Mutual_infos{i}),'r');
        plot(time,smooth(pos_Mutual_infos{i}),'k');
        
        xline(0)
        xlim([ -2300 1300])
        ylim([0 inf+0.1])
        xlabel('time shift')
        ylabel('MI')
        legend('Only blue','Only red','Original')
        if save_photo
           if i < 10
                 saveas(fig,[exp_folder,'\STA\FIG\',name,'_ch0',int2str(i),'.tif'])
            else
                saveas(fig,[exp_folder,'\STA\FIG\',name,'_ch',int2str(i),'.tif'])
           end
        end
        
%         figure('units','normalized','outerposition',[0 0 1 1])
%         ha = tight_subplot(1,2,[.04 .04],[0.09 0.02],[.04 .04]);
%         set(ha, 'Visible', 'off');
%         set(gcf,'units','normalized','outerposition',[0 0 1 1])
%         fig =gcf;
%         fig.PaperPositionMode = 'auto';
%         fig.InvertHardcopy = 'off';
%         axes(ha(1));
%         plot(time,smooth(PCA_Mutual_infos{1,i}-mean(PCA_Mutual_shuffle_infos{1,i})),'b');hold on
%         plot(time,smooth(PCA_Mutual_infos{2,i}-mean(PCA_Mutual_shuffle_infos{2,i})),'r');
%         plot(time,smooth(PCA_Mutual_infos{4,i}-mean(PCA_Mutual_shuffle_infos{4,i})),'b--');
%         plot(time,smooth(PCA_Mutual_infos{5,i}-mean(PCA_Mutual_shuffle_infos{5,i})),'r--');
%         legend('blue pos','red pos','blue v','red v')
%         xline(0)
%         xlim([ -2300 1300])
%         ylim([0 inf+0.1])
%         xlabel('time shift')
%         ylabel('MI')
%         
%         axes(ha(2));
%         [pd,x_values,y] = get_distribution(bin_pos');
%         [blue_pd,bluex_values,blue_y] = get_distribution(positive_before_pos{i});
%         [red_pd,redx_values,red_y] = get_distribution(negative_before_pos{i});
%         histogram(bin_pos,'Normalization','pdf');hold on
%         plot(x_values,y,'k')
%         plot(bluex_values,blue_y,'b')
%         plot(redx_values,red_y,'r')
%          if save_photo
%              if i < 10
%                  saveas(fig,[exp_folder,'\STA\FIG\MI_',name,'_ch0',int2str(i),'.tif'])
%              else
%                 saveas(fig,[exp_folder,'\STA\FIG\MI_',name,'_ch',int2str(i),'.tif'])
%              end
%         end
    end
end