% find reponse peaks of Gaussian pulse%% Gaussian pulse ON/OFF stimulus
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
BinningInterval = 0.02;
files=[56 50]; % select files


for z = 1:length(files)
%     clearvars -except SamplingRate files BinningInterval roi rr z all_file
    file = all_file(files(z)).name
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
   
            
    % stimulus formation
    t=1/SamplingRate:1/SamplingRate:length(a_data(1,:))/SamplingRate;
    sti=a_data(1,:);
    SingleSti_raw=sti(t>TimeStamps(1) & t<TimeStamps(2));
    tSingleSti=1/SamplingRate:1/SamplingRate:length(SingleSti_raw)/SamplingRate;
    % transform stimulus from raw value to volt
    SingleSti_volt=(SingleSti_raw-32768).*125*10^(-6); % transfrom from raw data of MCRack to output volt
    SingleSti=SingleSti_volt;
    SingleSti=smooth(SingleSti,100);
    % transform to real intensity
    load('F:\我的雲端硬碟\Retina exp\exp data\temporal_pattern\20200419\calibration_PAC_19-Apr-2020.mat')
    SingleSti=SingleSti-offset;
    Ip=SingleSti/10.421/10^6;
    r=0.37;
    P=Ip/r;
    A=13*10^-6;
    inten=P/A*1000; % unit: mW/m^2
    SingleSti=inten;
    
    
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
    
    % determine the stimulus peak
    t_stipk=2.5; % remember to change for different experiment 
%     meanSingleSti=SingleSti-mean(SingleSti);
%     [pkval stipklocs]=max(abs(meanSingleSti));
%     t_stipk=tSingleSti(stipklocs);
    rangedFR=cell(1,60);
    windowstart=t_stipk-0.5;
    windowend=t_stipk+0.5
    for i=1:60
        rangedFR{i}=FR{i}(tBinning>windowstart & tBinning<windowend);
    end
    rangedT=tBinning(tBinning>windowstart & tBinning<windowend);
    rangedTsti=tSingleSti(tSingleSti>windowstart & tSingleSti<windowend);
    rangedSti=SingleSti(tSingleSti>windowstart & tSingleSti<windowend);
%     figure
%     yyaxis left;plot(rangedT,rangedFR{32})
%     yyaxis right;plot(rangedTsti,rangedSti)
    
    % find apparent peak channel
    acceptchannel=[];
    for c=1:60
        peakrate=max(rangedFR{c})/mean(rangedFR{c});
        if peakrate>3
            if max(rangedFR{c})>10
                acceptchannel=[acceptchannel c];
            end
        end
    end
    
    load('F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\20200408_OU_original_tau=0p5_sort_unit2_MI_also_other_parameters.mat')
    P_accept=[];
    N_accept=[];
    for ii=1:length(acceptchannel)
        for jj=1:length(P_channel)
            if P_channel(jj)==acceptchannel(ii)
                P_accept=[P_accept acceptchannel(ii)];
                break
            end
        end
        for jj=1:length(N_channel)
            if N_channel(jj)==acceptchannel(ii)
                N_accept=[N_accept acceptchannel(ii)];
                break
            end
        end
    end
    
    %% find peak times
    % p channel
    for i=1:length(P_accept)
        channel=P_accept(i);
        [pkfrP(i) pklocP(i)]=max(rangedFR{channel});
    end
    % n channel
    for i=1:length(N_accept)
        channel=N_accept(i);
        [pkfrN(i) pklocN(i)]=max(rangedFR{channel});
    end
    
    %% multiplt array
%     colors={'k','r'};
    % P channel
    for n=1:length(P_accept)
        channel=P_accept(n);
        figure(z);subplot(8,8,rr(channel));hold on
        yyaxis left
        plot(rangedT,rangedFR{channel},'LineWidth',1,'LineStyle','-','Color','r');
        plot(rangedT(pklocP(n)),pkfrP(n),'ro')
        yyaxis right
        plot(rangedTsti,rangedSti,'LineWidth',1,'LineStyle','--','Color','k')
%         xlim([0 tSingleSti(end)-2])
    end
    
    % N channel
    for n=1:length(N_accept)
        channel=N_accept(n);
        figure(z);subplot(8,8,rr(channel));hold on
        yyaxis left
        plot(rangedT,rangedFR{channel},'LineWidth',1,'LineStyle','-','Color','b');
        plot(rangedT(pklocN(n)),pkfrN(n),'bo')
        yyaxis right
        plot(rangedTsti,rangedSti,'LineWidth',1,'LineStyle','--','Color','k')
%         xlim([0 tSingleSti(end)-2])
    end
    
    sign=['o','+'];
    figure(33);hold on
    plot(rangedT(pklocP),pkfrP,sign(z),'color','r');hold on
    plot(rangedT(pklocN),pkfrN,sign(z),'color','b')
    
    pkfrP_cell{z}=pkfrP;
    pkfrN_cell{z}=pkfrN;
end
