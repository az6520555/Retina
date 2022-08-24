% Find prediction horizon under different sitmulation intensity
clear
close all
    rr =[9,17,25,33,41,49,...
      2,10,18,26,34,42,50,58,...
      3,11,19,27,35,43,51,59,...
      4,12,20,28,36,44,52,60,...
      5,13,21,29,37,45,53,61,...
      6,14,22,30,38,46,54,62,...
      7,15,23,31,39,47,55,63,...
        16,24,32,40,48,56];
path=['F:\我的雲端硬碟\Retina exp\exp data\Sorted_final_data\20200408\MIdata'];
load(['F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\20200408_OU_original_tau=0p5_sort_unit2_MI_also_other_parameters.mat']); % import on off classification
P_channel=[1,20,23,33,39,40,42,53];
P_channel_all=P_channel; % the parameter name is confront with another loaded files 
N_channel_all=N_channel; % the parameter name is confront with another loaded files 
cd(path);
all_file = dir('*.mat');
file_seq=[2 11 8 5]; 
    % 9 7 5 13 for 20210331
    % 6 8 9 13 for 20210413
    % 5 9 11 18 13 15 7 for 20210420\
barinten_seq=[5,3.5,2,1]; 
    % 5.4 4.333 3.6111 0 for 20210331 (background intensity and 6.5 bar intensity)
    % 1.3 3.25 5.2 6.5 for 20210413 (bar intensity with zero background intensity)
    % 1.3 3.25 5.2 6.5 7.8 9.75 11.7 for 20210420 (bar intenisty)
roi=1:60;
% ====subplot axis titles======
x_title='cut-off frequency (Hz)'; %'Bar intensity (mW/m^2)' ,'Background intensity (mW/m^2)'
y_title='\deltat_{peak} (ms)';

cc={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]};

for z = 1:length(file_seq)
    load(fullfile(path,all_file(file_seq(z)).name))
%     TimeShift=time;
%     MI=Mutual_infos;
    % exclude MI
    for ch=roi
        [MI_peak(z,ch) ind_peak(z,ch)] = max(MI{ch}); % MI_peak{number of data}{number of channel}
        if (TimeShift(ind_peak(z,ch)) < -500) || (TimeShift(ind_peak(z,ch)) > 1000) % the 100 points, timeshift is 1000
            MI_peak(z,ch) = NaN;
            ind_peak(z,ch) = NaN;
        end
    end
    % exclude the MI that is too small 
    for nn=roi
        if MI_peak(z,nn)-min(MI{nn}) < 0.15
            MI_peak(z,nn) = NaN;
            ind_peak(z,nn) = NaN;
            MI{nn}=[];
        end
    end 
end
%% calculate MI peak 
for i=1:size(ind_peak,1)
    for j=1:size(ind_peak,2)
        if isnan(ind_peak(i,j))
            all_MI_peakshift(i,j)=0;
            all_MI_peakshift2(i,j)=NaN;
            ispk(i,j)=0;
        else
            all_MI_peakshift(i,j)=TimeShift(ind_peak(i,j));
            all_MI_peakshift2(i,j)=all_MI_peakshift(i,j);
            ispk(i,j)=1;
        end
    end
end

% ======= mean error of N, P cells
npks_P=sum(ispk(:,P_channel_all),2); % 1*n array, number of cells which are detected MI peaks under n different bar intensity
npks_N=sum(ispk(:,N_channel_all),2); % 1*n array, number of cells which are detected MI peaks under n different bar intensity
for i=1:size(ind_peak,1)
    P_mean_dt(i)=sum(all_MI_peakshift(i,P_channel_all))/npks_P(i); 
    P_var_dt(i)=std(all_MI_peakshift2(i,P_channel_all),'omitnan');
    N_mean_dt(i)=sum(all_MI_peakshift(i,N_channel_all))/npks_N(i);
    N_var_dt(i)=std(all_MI_peakshift2(i,N_channel_all),'omitnan');
end
figure(1);sgtitle('MI peak timelag of N and P cells')
subplot(2,2,1) % plot MI peak position separately
for p_ch=P_channel_all
    try
        
        plot(barinten_seq,TimeShift(ind_peak(:,p_ch)),'-o','linewidth',1);hold on
    catch
    end
end 

title('Single P type cells')
xlabel(x_title)
ylabel(y_title)
subplot(2,2,2) % mean peak position
errorbar(barinten_seq,P_mean_dt,P_var_dt,'o','color','r','linewidth',1)
title('Mean and error of P type cells')
xlabel(x_title)
ylabel(y_title)
subplot(2,2,3) % plot peak position separately
for n_ch=N_channel_all
    try
        plot(barinten_seq,TimeShift(ind_peak(:,n_ch)),'-o','linewidth',1);hold on
    catch
    end
end 

title('Single N type cells')
xlabel(x_title)
ylabel(y_title)
subplot(2,2,4)
errorbar(barinten_seq,N_mean_dt,N_var_dt,'o','color','k','linewidth',1)
title('Mean and error of N type cells')
xlabel(x_title)
ylabel(y_title)

%% set first timelag as zero
MI_peakshift_1zero=all_MI_peakshift-repmat(all_MI_peakshift(1,:),size(all_MI_peakshift,1),1);
figure(2);sgtitle('MI peak timelag of N and P cells (subtract timelag of the lowest contrast)')
subplot(2,2,1) % P type
for ch=P_channel_all
    if all(ispk(:,ch))
        plot(barinten_seq,MI_peakshift_1zero(:,ch),'-o','linewidth',1);hold on
    end
end
yline(0,'--')
title('Single P type cells')
xlabel(x_title)
ylabel(y_title)
subplot(2,2,3)
for ch=N_channel_all % N type
    plot(barinten_seq,MI_peakshift_1zero(:,ch),'-o','linewidth',1);hold on
end
yline(0,'--')
title('Single N type cells')
xlabel(x_title)
ylabel(y_title)

for i=1:size(ind_peak,1)
    P_mean_dt_1zero(i)=sum(MI_peakshift_1zero(i,P_channel_all))/npks_P(i);
    P_var_dt_1zero(i)=std(MI_peakshift_1zero(i,P_channel_all));
    N_mean_dt_1zero(i)=sum(MI_peakshift_1zero(i,N_channel_all))/npks_N(i);
    N_var_dt_1zero(i)=std(MI_peakshift_1zero(i,N_channel_all));
end
subplot(2,2,2)
errorbar(barinten_seq,P_mean_dt_1zero,P_var_dt_1zero,'o','color','r','linewidth',1);
title('Mean and error of P type cells')
xlabel(x_title)
ylabel(y_title)
subplot(2,2,4)
errorbar(barinten_seq,N_mean_dt_1zero,N_var_dt_1zero,'o','color','k','linewidth',1);
title('Mean and error of N type cells')
xlabel(x_title)
ylabel(y_title)

%% plot only errorbar
figure(111);hold on;box on;
errorbar(barinten_seq,P_mean_dt,P_var_dt,'o','color','r','linewidth',1)
errorbar(barinten_seq,N_mean_dt,N_var_dt,'o','color','k','linewidth',1)
% title('Mean and error of N and NP type cells')
xlabel(x_title)
ylabel(y_title)
yline(0,'--','color',[0.5,0.5,0.5])
legend('P cell','NP cell')
xlim([0.5,5.5])