%% Gaussian pulse ON/OFF stimulus
clear all
close all

path=['\\192.168.0.102\Public\Retina\Chou\Exp\data_until_2020\Sorted_final_data\20200408'];
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
files=[1 10]; % 25 31 28 34

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
   

    %% stimulus formation
    t=1/SamplingRate:1/SamplingRate:length(a_data(1,:))/SamplingRate;
    sti=a_data(1,:);
    SingleSti_raw=sti(t>TimeStamps(1) & t<TimeStamps(2));
    tSingleSti=1/SamplingRate:1/SamplingRate:length(SingleSti_raw)/SamplingRate;
    %% transform stimulus from raw value to volt
    SingleSti_volt=(SingleSti_raw-32768).*125*10^(-6); % transfrom from raw data of MCRack to output volt
    SingleSti=SingleSti_volt;
    SingleSti=smooth(SingleSti,100);
    %% transform to real intensity
%     load('F:\我的雲端硬碟\Retina exp\exp data\temporal_pattern\20200419\calibration_PAC_19-Apr-2020.mat')
%     SingleSti=SingleSti-offset;
%     Ip=SingleSti/10.421/10^6;
%     r=0.37;
%     P=Ip/r;
%     A=13*10^-6;
%     inten=P/A*1000; % unit: mW/m^2
%     SingleSti=inten;
    
    
    %% Binning spikes
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
    
    %% find gaussian peak time and firing peak time
    [M,I]=max(abs(SingleSti-10));
    gaussian_peaktime=tSingleSti(I);
    
    for i=1:60
        [FR_peak(i) FR_peakind(i)]=max(FR{i});
        FR_peaktime(i)=tBinning(FR_peakind(i));
    end
    
%%     ========plot multiple arrays (add on off index)============
    colors={'k','r'};
    figure(1)
    for n=1:60
        subplot(8,8,rr(n));hold on
        yyaxis left
        plot(tBinning,FR{n},'LineWidth',1,'LineStyle','-','Color',colors{z});
        yyaxis right
        plot(tSingleSti,SingleSti,'LineWidth',1,'LineStyle','--','Color',colors{z})
%         xlim([0 tSingleSti(end)-2])
        xlim([0 2])
    end
%     pathonoff='F:\我的雲端硬碟\Retina exp\exp data\整理\OnOff_index\';
%     file_onoff=[pathonoff,'OU_tau=600ms_19-Apr-2020_0_sort_unit1_MI_also_other_parameters.mat'];
%     onoff_color(file_onoff)
%     pathNP='F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\';
%     file_NP=[pathNP,'Gonoff-19-Apr-2020-0-sort-unit1onoff_index.mat'];
%     NPcolor(file_NP)

    %% =========== All channel raster plot===========
%     for i=1:60
%         if isempty(Spikes{i})==1
%             Spikes{i}=0;
%         end
%     end
%     figure(z*2);
%     LineFormat.Color = [0 0 0];
%     subplot(5,1,[1:4])
%     plotSpikeRaster(Spikes,'PlotType','vertline','LineFormat',LineFormat)
%     subplot(5,1,5)
%     plot(t,sti)

%     samexaxis('abc','xmt','on','ytac','join','yld',1);

    %% =============compare channel of on and off=============
    figure;hold on
    channel=39;
    title(['channel ',num2str(channel),'  response delay = ',num2str(FR_peaktime(channel)-gaussian_peaktime),' second'])
%     subplot(length(files),1,z)
    yyaxis left
    plot(tBinning,FR{channel},'LineWidth',2);%'LineStyle','-'
    ylabel('Firing Rate (Hz)')
%     ylim([0 60])
    yyaxis right
    plot(tSingleSti,SingleSti,'LineWidth',2)
    ylabel('Stimuli (mW/m^2)')
%     xlim([0 tSingleSti(end)-2])
    xlim([2,4])
    disp(['response delay = ',num2str(FR_peaktime(channel)-gaussian_peaktime),' second'])
    
    
%     fr_save=[fr_save tBinning(tBinning>1 & tBinning<5)'-1 FR{channel}(tBinning>1 & tBinning<5)'];
%     sti_save=[sti_save tSingleSti(tSingleSti>1 & tSingleSti<5)'-1 SingleSti(tSingleSti>1 & tSingleSti<5)];
    
%% classify the behavior of on, on-off, off cell
%     load('F:\我的雲端硬碟\Retina exp\exp data\整理\OnOff_index\Gonoff-20-Apr-2020-0-sort-unit1onoff_index.mat')
%     
%     colorsNF={'k','r'};
%     figure(333);hold on;box on
%     %ON
%     subplot(7,1,[1,2]);hold on;grid on;box on
%     onsum=zeros(1,length(FR{1}));
%     for i=1:length(On_channel)
%         onsum=onsum+FR{On_channel(i)};
%     end
%     onsum=onsum/length(On_channel);
%     plot(tBinning,onsum,colorsNF{z})
%     ylabel({['ON cell  N=',num2str(length(On_channel))];'Firing Rate (hz)'})
%     
%     %OFF
%     offsum=zeros(1,length(FR{1}));
%     subplot(7,1,[3,4]);hold on;grid on;box on
%     for i=1:length(Off_channel)
%         offsum=offsum+FR{Off_channel(i)};
%     end
%     offsum=offsum/length(Off_channel);
%     plot(tBinning,offsum,colorsNF{z})
%     ylabel({['OFF cell  N=',num2str(length(Off_channel))];'Firing Rate (hz)'})
%     
%     %ONOFF
%     onoffsum=zeros(1,length(FR{1}));
%     subplot(7,1,[5,6]);hold on;grid on;box on
%     for i=1:length(OnOff_channel)
%         onoffsum=onoffsum+FR{OnOff_channel(i)};
%     end
%     onoffsum=onoffsum/length(OnOff_channel);
%     plot(tBinning,onoffsum,colorsNF{z})
%     ylabel(['ON-OFF  N=',num2str(length(OnOff_channel))])
%     
%     subplot(7,1,7);hold on;grid on;box on
%     plot(tSingleSti,SingleSti,'LineWidth',1.5,'color',colorsNF{z})
%     ylabel('Stimulus (mW/m^{2})')
%     xlabel('time (s)')
%     
    %% classify the behavior of predictive and non-predictive cell
%     load('F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\OU_tau=100ms_19-Apr-2020_0_sort_unit1_MI_also_other_parameters.mat')
%     colorsNP={'b','g'};
%     figure(877);hold on
%     
%     subplot(5,1,[1,2]);hold on;grid on;box on
%     Nsum=zeros(1,length(FR{1}));
%     for i=1:length(N_channel)
%         Nsum=Nsum+FR{N_channel(i)};
%     end
%     Nsum=Nsum/length(N_channel);
%     plot(tBinning,Nsum,colorsNP{z})
%     ylabel({['N type cell  N=',num2str(length(N_channel))];'Firing Rate (hz)'})
%     
%     subplot(5,1,[3,4]);hold on;grid on;box on
%     Psum=zeros(1,length(FR{1}));
%     for i=1:length(P_channel)
%         Psum=Psum+FR{P_channel(i)};
%     end
%     Psum=Psum/length(P_channel);
%     plot(tBinning,Psum,colorsNP{z})    
%     ylabel({['P type cell  N=',num2str(length(P_channel))];'Firing Rate (hz)'})
%     
%     subplot(5,1,5);hold on;grid on;box on
%     plot(tSingleSti,SingleSti,'LineWidth',1.5,'color',colorsNP{z})
%     ylabel('Stimulus (mW/m^{2})')
%     xlabel('time (s)')
%     
%     
%     figure(222);hold on;grid on
%     plot(tSingleSti,SingleSti,'LineWidth',1.5,'color',colorsNF{z})
%     legend('ON','OFF')
%     ylabel('Stimulus (mW/m^{2})')
%     xlabel('time (s)')

%% output mat data
% matname={'gaussianON.mat','gaussianOFF.mat'};
% frN=FR{32};
% save(['F:\我的雲端硬碟\Retina exp\Retina as NGD paper submit\figure 20201019 modification\',matname{z}] ...
%     ,'tBinning','frN','tSingleSti','SingleSti')
end
% figure(287);samexaxis('abc','xmt','on','ytac','join','yld',1);
% try
%     figure(333);samexaxis('abc','xmt','on','ytac','join','yld',1);
% catch
% end
% try
%     figure(877);samexaxis('abc','xmt','on','ytac','join','yld',1);
% catch
% end
% 
%% save paper data
% sti_save=sti_save(2:2:end,:); % prevent the matrix exceeding the limit
% path_data='F:\我的雲端硬碟\Retina exp\Retina as NGD paper submit\data';
% data_name_sti='gaussian_pw=0p6_stimulus.xls';
% writematrix(sti_save,fullfile(path_data,data_name_sti))
% data_name_fr='gaussian_pw=0p6_response.xls';
% writematrix(fr_save,fullfile(path_data,data_name_fr))




