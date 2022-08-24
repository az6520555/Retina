% recognize N type or P type cell
% load MI data
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
Date='20210513';
path=['F:\我的雲端硬碟\Retina exp\exp data\Spatial stimuli\',Date,'\MI\unsort'];
% load('F:\我的雲端硬碟\Retina exp\exp data\整理\OnOff_index\Gonoff-19-Apr-2020-0-sort-unit1onoff_index'); % import on off classification
cd(path);
all_file = dir('*.mat');
file_seq=[20 19 18 17];
% corr_time=[[0.12285,0.3345,0.43395,0.7175,1.27435]];
roi=1:60;
cc={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]};

for z = 1:length(file_seq)
    load(fullfile(path,all_file(file_seq(z)).name))
    TimeShift=time;
    MI=Mutual_infos;
    for ch=roi
        [MI_peak(z,ch) ind_peak(z,ch)] = max(MI{ch}); % MI_peak{number of data}{number of channel}
        if (TimeShift(ind_peak(z,ch)) < -500) || (TimeShift(ind_peak(z,ch)) > 1000) % the 100 points, timeshift is 1000
            MI_peak(z,ch) = NaN;
            ind_peak(z,ch) = NaN;
            MI{ch}=[];
        end
    end
    % ======= exclude the MI that is too small ===========
    for nn=roi
        if MI_peak(z,nn)-min(MI{nn}) < 0.15 
            MI_peak(z,nn) = NaN;
            ind_peak(z,nn) = NaN;
            MI{nn}=[];
        end
    end
    %% plot multichannel MI
    figure(1);
    if z==1
        for i=1:60
            subplot(8,8,rr(i));hold on
            xline(0,':') 
        end
    end
    for nn=roi
        subplot(8,8,rr(nn));hold on
        title(num2str(nn))
        try
            plot(TimeShift,MI{nn},'LineWidth',2,'LineStyle','-','color',cc{z},'DisplayName',all_file(file_seq(z)).name(1:end-4));
            plot(TimeShift(ind_peak(z,nn)),MI_peak(z,nn),'Marker','d','color',cc{z},'DisplayName',all_file(file_seq(z)).name(1:end-4))
        catch
            continue
        end
        xlim([-1000 1000])
    end
    legend('Interpreter', 'none', 'Position', [0.89 0.92 0.03 0.02])
    
    %% plot peak shift
%     for nn=roi
%         figure(2);subplot(8,8,rr(nn));hold on
%         
%         try
%             plot(corr_time(z),TimeShift(ind_peak(z,nn)),'-o','color',cc{z})
%         catch
%             continue
%         end
%         xlim([0 1.3])
%         ylim([-1000 1000])
%     end
    
    %% single plot
    channel=46;
    figure(5);hold on
    title('OU','FontSize',20)
    try 
        plot(TimeShift,MI{channel},'LineWidth',2,'LineStyle','-');
    catch
    end
    ylabel('MI[\gamma(t),I(t-\deltat)]','FontSize',18)
    xlabel('time shift (ms)','FontSize',18)
    refline([0 0])
    xlim([-1000 1000])
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'off';
end



%% separate the cells of different response (classify two types of cell)
% consider different number of peak point
% ======= exclude the channel that all MI are small (<0.15 bits/s)======


P_channel=[]; % predictive cell
N_channel=[]; % none predictive cell
Pslight_channel=[]; % slightly predictive cell
for n = roi
    pks_1d=ind_peak(:,n);
    peaks_ind=pks_1d(~isnan(pks_1d));
    peaks_time = TimeShift(peaks_ind); 
    % ============ find predictive cell=============
    pk_pos=peaks_time(peaks_time>0);
    if length(pk_pos)>=1
        P_channel=[P_channel n];
    end

    % ============ find none predictive cell ===============
    if (length(peaks_time)>=2) & (peaks_time<0)
        N_channel=[N_channel n];
    end
end

[pathstr, name, ext] = fileparts(all_file(file_seq(1)).name);
path_NP='F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\';
cd([path_NP,'NP_spatial'])
save([Date,'_NP','.mat'],'P_channel','N_channel','Pslight_channel')

%% change figure color to classify predictive or non-predictive cells
for i=[1 2]
    figure(i)
    % === P ===
    for n=1:length(P_channel)
        subplot(8,8,rr(P_channel(n)))
        set(gca,'Color',[0.8 1 0.8])
    end
    % === N ===
    for n=1:length(N_channel)
        subplot(8,8,rr(N_channel(n)))
        set(gca,'Color',[0.8 0.8 1])
    end
end

% ========== save figures ============
cd(path)
mkdir MIfigure_NP
figure(1)
set(gcf, 'InvertHardCopy', 'off');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
tic
n=0;
while n<0.5
    n=toc;
end
% figure(1)
% saveas(gcf,[path,'\MIfigure_NP\',Date,'_NP.fig']);
% print(gcf,[path,'\MIfigure_NP\',Date,'_NP.png'],'-dpng','-r300');

%% change axes color to classify on, on-off or off cells
% for ii=[1,2]
%     figure(ii)
%     for jj=1:length(On_channel)
%         subplot(8,8,rr(On_channel(jj)));hold on
%         box on
%         ax=gca;
%         ax.XColor=[0.6 0 0]; % red
%         ax.YColor=[0.6 0 0]; % red
%         ax.LineWidth = 1;
%     end
%     for jj=1:length(Off_channel)
%         subplot(8,8,rr(Off_channel(jj)));hold on
%         box on
%         ax=gca;
%         ax.XColor=[0 0.6 0]; % green
%         ax.YColor=[0 0.6 0]; % green
%         ax.LineWidth = 1;
%     end
%     for jj=1:length(OnOff_channel)
%         subplot(8,8,rr(OnOff_channel(jj)));hold on
%         box on
%         ax=gca;
%         ax.XColor=[0.6 0 0.6]; % magenta
%         ax.YColor=[0.6 0 0.6]; % magenta
%         ax.LineWidth = 1;
%     end
% end

%% calculate the mean shift of different type of cells

% for i=1:size(ind_peak,1)
%     nnn=0;
%     for j=1:size(ind_peak,2)
%         if isnan(ind_peak(i,j))==1
%             new_dt(i,j)=0;
%             ispk(i,j)=0;
%         else
%             new_dt(i,j)=TimeShift(ind_peak(i,j));
%             ispk(i,j)=1;
%         end
%     end
% end
% 
% % ====== for P type cell ===========
% npks_P=sum(ispk(:,P_channel),2);
% for i=1:size(ind_peak,1)    
%     P_mean_dt(i)=sum(new_dt(i,P_channel))/npks_P(i);
% end
% for i=1:size(ind_peak,1)
%     P_var_dt(i)=std(new_dt(i,P_channel));
% end
% figure(3)
% errorbar(corr_time,P_mean_dt,P_var_dt,'o','color','r');hold on;
% title('HMM Mutual information time shift of P type cells')
% xlabel('Correlation time (ms)')
% ylabel('Timeshift of MI peaks')
% 
% % ====== for N type cell ===========
% npks_N=sum(ispk(:,N_channel),2);
% for i=1:size(ind_peak,1)    
%     N_mean_dt(i)=sum(new_dt(i,N_channel))/npks_N(i);
% end
% for i=1:size(ind_peak,1)
%     N_var_dt(i)=std(new_dt(i,N_channel));
% end
% figure(4)
% errorbar(corr_time,N_mean_dt,N_var_dt,'o','color','k')
% title('HMM Mutual information time shift of N type cells')
% xlabel('Correlation time (ms)')
% ylabel('Timeshift of MI peaks')


% figure(3);saveas(gcf,[path,'MIfig\',name,'_and_other_different_parameters_Ppeaks.png']);
% figure(4);saveas(gcf,[path,'MIfig\',name,'_and_other_different_parameters_Npeaks.png']);


