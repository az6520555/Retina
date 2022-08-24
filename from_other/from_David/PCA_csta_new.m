%% This code calculate HMM and OU bar position STA
close all;
clear all;
exp_folder = 'F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200419\';
cd(exp_folder);
save_photo =1;%0 is no save on off photo and data, 1 is save
file = 'OU_tau=600ms_cutoff=4_19-Apr-2020_0_sort_unit1.mat';%Name that used to save photo and data 'OU_tau=100ms_cutoff=2_08-May-2020_0_sort_unit1.mat'
[pathstr, name, ext] = fileparts(file);
bin = 10;  BinningInterval = bin*10^-3;  %ms
mkdir PCAdata
% %% For unsorted spikes
% load(['data\',name,'.mat'])
%% For sorted spikes
load([exp_folder,file])
SamplingRate = 20000;
%% a_data as TriggerData
TriggerData = a_data(1,TimeStamps(1)*SamplingRate:TimeStamps(length(TimeStamps))*SamplingRate);% figure;plot(isi);
inten = downsample(TriggerData,SamplingRate*BinningInterval);
inten = inten-mean(inten);
% inten = 2*(inten-min(inten))/(max(inten)-min(inten))-1;
%%


sorted = 1;
unit = 1;
forward = 90;%90 bins before spikes for calculating STA
backward = 90;%90 bins after spikes for calculating STA
%roi = [p_channel np_channel];
roi = 39;

%     try
%         load([exp_folder,'\STA\',name(12:end),'.mat'])
%     catch
idxs = cell(1,60);
PCA_STAs = cell(1,60);
positive_PCAs = cell(1,60);
negative_PCAs = cell(1,60);
BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
BinningSpike = zeros(60,length(BinningTime));
positive_before_pos = cell(1,60);
negative_before_pos = cell(1,60);
%     end
TheStimuli= inten;  %recalculated bar position
acf = autocorr(TheStimuli,100);
corr_time = interp1(acf,1:length(acf),0.5,'linear')/60;
   %%  Binning
bin=BinningInterval*10^3; %ms

STA_time=[-forward:backward]*BinningInterval;
%% BinningSpike and calculate STA
analyze_spikes = Spikes;
sum_n = zeros(1,60);
dis_STA = zeros(60,forward+backward+1);
mkdir([exp_folder,'PCAdata\PCAfig\',name]) % directory to save PCA diiagram
for i = roi  % i is the channel number
    if length(analyze_spikes{i}) > 5
        [n,~] = hist(analyze_spikes{i},BinningTime);
        PCA_STA = zeros(1,forward+backward+1);
        BinningSpike(i,:) = n ;
        sum_n(i) = sum_n(i)+ sum(BinningSpike(i,forward+1:length(BinningTime)-backward));
        num_spike =1;
        pos_100ms = [];
        for time_shift = forward+1:length(BinningTime)-backward
            dis_STA(i,:) = dis_STA(i,:) + BinningSpike(i,time_shift)*TheStimuli(time_shift-forward:time_shift+backward);
            for x = 1:BinningSpike(i,time_shift)
                pos_100ms = [pos_100ms;TheStimuli(time_shift-6)];
                PCA_STA(num_spike,:) = TheStimuli(time_shift-forward:time_shift+backward);
                num_spike = num_spike+1;
            end

        end
        if sum_n(i)
            dis_STA(i,:) = dis_STA(i,:)/sum_n(i);
        end
        [coeff,score,latent,tsquared,explained] = pca(PCA_STA);
        PCA_STAs{i} = PCA_STA;
        %% K means clusterin2
        score1 = score(:,1);
        score2 = score(:,2);
        idx = zeros(1,length(score));
        idx(find(score1<0)) = 1;
        idx(find(score1>=0)) = 2;
%         cluster_data = [score1,score2];%Merge PCA1 and PCA2
%         idx = kmeans(cluster_data,2);%Results of which clusters(1 or 2), and clustered into 2 groups
        idxs{i} = idx;
        
    %         %Red represents group 1, blue represents group 2
        figure(i)
        scatter(score1(find(idx==1)),score2(find(idx==1)),'b');hold on
        scatter(score1(find(idx==2)),score2(find(idx==2)),'r');
        xlabel('Stimulus*PCA1')
        ylabel('Stimulus*PCA2')
        saveas(gcf,[exp_folder,'PCAdata\PCAfig\',name,'\PCA_channel',num2str(i),'.png'])
        title('Results of two cluster')
        
    %     axes(ha(2)); 
        positive_PCA = PCA_STA(find(idx==1),:);
        positive_PCAs{i} = positive_PCA;
        negative_PCA= PCA_STA(find(idx==2),:);
        negative_PCAs{i} = negative_PCA;
        positive_before_pos{i} = pos_100ms(find(idx==1));
        negative_before_pos{i} = pos_100ms(find(idx==2));
        
% 
        figure;hold on
        STA_sum=mean(PCA_STA,1);
        STA1=mean(positive_PCAs{i},1);
        STA2=mean(negative_PCAs{i},1);
        plot(STA_time,STA_sum,'k')
        plot(STA_time,STA1,'b')%Red is positive_PCA
        plot(STA_time,STA2,'r');%Blue is negative_PCA
        yline(0)
        legend('STA','STA cluster 1','STA cluster 2')
        saveas(gcf,[exp_folder,'PCAdata\PCAfig\',name,'\STA_channel',num2str(i),'.png'])
        
%         figure;plot(diff(STA2)*2);hold on;plot(STA1);plot(STA2)
%         figure;hold on;plot(STA_time,STA1,'b');plot(STA_time,-STA2,'r')
    end
end
save([exp_folder,'\PCAdata\',file],'STA_time','dis_STA','corr_time','idxs','PCA_STAs','positive_PCAs',...
'negative_PCAs','positive_before_pos','negative_before_pos','BinningSpike','TheStimuli','BinningTime','forward','backward','bin','BinningInterval','score1','score2')


