%% Gaussian pulse ON/OFF stimulus
clear all
close all

path=['F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200419'];
cd(path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
SamplingRate=20000;
cc=hsv(n_file);

        rr =[9,17,25,33,41,49,...
          2,10,18,26,34,42,50,58,...
          3,11,19,27,35,43,51,59,...
          4,12,20,28,36,44,52,60,...
          5,13,21,29,37,45,53,61,...
          6,14,22,30,38,46,54,62,...
          7,15,23,31,39,47,55,63,...
            16,24,32,40,48,56];
roi = 1:60;
BinningInterval = 0.01;
files=[]; % 34 28 25 31

sti_save=[];
fr_save=[];


for z = 1:length(files)
    clearvars -except SamplingRate files BinningInterval roi rr z all_file sti_save fr_save
    file = all_file(files(z)).name 
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
   

    % stimulus formation
    t=1/SamplingRate:1/SamplingRate:length(a_data(1,:))/SamplingRate;
    sti=a_data(1,:);
    SingleSti_raw=sti(t>TimeStamps(5) & t<TimeStamps(6));
    tSingleSti=1/SamplingRate:1/SamplingRate:length(SingleSti_raw)/SamplingRate;
    % transform stimulus from raw value to volt
    SingleSti_volt=(SingleSti_raw-32768).*125*10^(-6); % transfrom from raw data of MCRack to output volt
    SingleSti=SingleSti_volt;
    SingleSti=smooth(SingleSti,1000);
    % transform to real intensity
    load('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\temporal_pattern\20200419\calibration_PAC_19-Apr-2020.mat')
    SingleSti=SingleSti-offset;
    Ip=SingleSti/10.421/10^6;
    r=0.37;
    P=Ip/r;
    A=13*10^-6;
    inten=P/A*1000; % unit: mW/m^2
    SingleSti=inten;
    c=['k','r']
    
%     figure(8787);hold on;box on
%     plot(tSingleSti,SingleSti,'linewidth',1,'color',c(z))
%     title('Gaussian pulse')
%     xlim([1 5])
%     ylim([4 16])
%     xlabel('time (s)')
%     ylabel('light intensity (mW/m^2)')
%     legend('ON','OFF')
%     set(gcf,'Position',[300,300,300,200])
    
    
    % Binning spikes
    trial=length(TimeStamps)-1;
    SpikesSupPos=cell(1,60);
    for i=1:60
        for j=1:trial-1
            sp=[];
            SpikesST=[];
            sp=Spikes{i};
            SpikesST=sp(sp>TimeStamps(j) & sp<TimeStamps(j+1))-TimeStamps(j);
            SpikesSupPos{i}=cat(2,SpikesSupPos{i},SpikesST);
        end
    end
    BinningTime=0:BinningInterval:TimeStamps(2)-TimeStamps(1);
    FR=cell(1,60);
    for i=1:60
        [counts,centers]=hist(SpikesSupPos{i},BinningTime);
        FR{i}=counts/(trial-1)/BinningInterval;
    end
    tBinning=centers;
    
    %% plot three channels
    
    %% find gaussian peak time and firing peak time
%     [M,I]=max(abs(SingleSti-10));
%     gaussian_peaktime=tSingleSti(I);
%     
%     for i=1:60
%         [FR_peak(i) FR_peakind(i)]=max(FR{i});
%         FR_peaktime(i)=tBinning(FR_peakind(i));
%     end
    


    %% =============compare channel of on and off=============
%     figure;hold on
%     channel=39;
%     title(['channel ',num2str(channel),'  response delay = ',num2str(FR_peaktime(channel)-gaussian_peaktime),' second'])
% %     subplot(length(files),1,z)
%     yyaxis left
%     plot(tBinning,FR{channel},'LineWidth',2);%'LineStyle','-'
%     ylabel('Firing Rate (Hz)')
% %     ylim([0 60])
%     yyaxis right
%     plot(tSingleSti,SingleSti,'LineWidth',2)
%     ylabel('Stimuli (mW/m^2)')
% %     xlim([0 tSingleSti(end)-2])
%     xlim([2,4])
%     disp(['response delay = ',num2str(FR_peaktime(channel)-gaussian_peaktime),' second'])
%     
%     
%     fr_save=[fr_save tBinning(tBinning>1 & tBinning<5)'-1 FR{channel}(tBinning>1 & tBinning<5)'];
%     sti_save=[sti_save tSingleSti(tSingleSti>1 & tSingleSti<5)'-1 SingleSti(tSingleSti>1 & tSingleSti<5)];
%     


%% output mat data
% frN={FR{46},FR{10},FR{25}};
% save(['F:\§Úªº¶³ºÝµwºÐ\Master Thesis\Figures\Result\data\',name,'.mat'] ...
%     ,'tBinning','frN','tSingleSti','SingleSti')
end



