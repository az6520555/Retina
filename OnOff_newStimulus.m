%% analysis onoff (newone)
clear all
close all

path='F:\我的雲端硬碟\Retina exp\exp data\Sorted_final_data\20200318';
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


for z = []
    file = all_file(z).name;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    % remove low firing rate
    for i=1:60
        if length(Spikes{i})/(TimeStamps(end)-TimeStamps(1))<1
            Spikes{i}=[];
        end
    end
            
    % stimulus formation
    t=1/SamplingRate:1/SamplingRate:length(a_data(1,:))/SamplingRate;
    sti=a_data(1,:);
    SingleSti1=sti(t>TimeStamps(1) & t<TimeStamps(1)+1);
    SingleSti2=sti(t>TimeStamps(2) & t<TimeStamps(2)+1);
    SingleSti1=SingleSti1-min(SingleSti1); 
    SingleSti2=SingleSti2-min(SingleSti2); 
    tSingleSti1=1/SamplingRate:1/SamplingRate:length(SingleSti1)/SamplingRate;
    tSingleSti2=1/SamplingRate:1/SamplingRate:length(SingleSti2)/SamplingRate;
    
    % Binning spikes
    trial=length(TimeStamps)/2;
    SpikesSupPos1=cell(1,60);
    SpikesSupPos2=cell(1,60);
    for i=1:60
        for j=1:trial
            sp=[];
            SpikesST1=[];
            SpikesST2=[];
            sp=Spikes{i};
            SpikesST1=sp(sp>TimeStamps(2*j-1) & sp<TimeStamps(2*j-1)+1)-TimeStamps(2*j-1);
            SpikesST2=sp(sp>TimeStamps(2*j) & sp<TimeStamps(2*j)+1)-TimeStamps(2*j);
            SpikesSupPos1{i}=cat(2,SpikesSupPos1{i},SpikesST1);
            SpikesSupPos2{i}=cat(2,SpikesSupPos2{i},SpikesST2);
        end
    end
    
    %===========calculate On-Off index=================
    onoff_index=[];
    for ii=1:60
        s1=[];
        s2=[];
        s1=SpikesSupPos1{ii};
        s1=s1(s1>0.05 & s1<0.55);
        s2=SpikesSupPos2{ii};
        s2=s2(s2>0.05 & s2<0.55);
        onoff_index(ii)=(length(s1)-length(s2))/(length(s1)+length(s2));
    end
    % ========classifying on, on-off or off cell============
    On_channel=[];
    Off_channel=[];
    OnOff_channel=[];
    for k=1:60
        if onoff_index(k)>0.3
            On_channel=[On_channel k];
        elseif (onoff_index(k)<=0.3) && (onoff_index(k)>=-0.3)
            OnOff_channel=[OnOff_channel k];
        elseif onoff_index(k)<-0.3
            Off_channel=[Off_channel k];
        end
    end
    % save on off index
    onoff_path='F:\我的雲端硬碟\Retina exp\exp data\整理\OnOff_index\';
    save([onoff_path,name,'onoff_index.mat'],'On_channel','OnOff_channel','Off_channel','onoff_index')
    %==========================================
    
    BinningTime=0:BinningInterval:1;
    FR1=cell(1,60);
    for i=1:60
        [counts1,centers]=hist(SpikesSupPos1{i},BinningTime);
        [counts2,centers]=hist(SpikesSupPos2{i},BinningTime);
        FR1{i}=counts1/(trial)/BinningInterval;
        FR2{i}=counts2/(trial)/BinningInterval;
    end
    tBinning=centers;
    
%     %============ plot multiple arrays==================
%     for n=1:60
%         figure(z);subplot(8,8,rr(n));hold on
%         yyaxis left
%         plot(tBinning,FR1{n},'LineWidth',2,'LineStyle','-');
%         yyaxis right
%         plot(tSingleSti1,SingleSti1,'LineWidth',2,'LineStyle','-')
%     end
%     for n=1:60
%         figure(z*2);subplot(8,8,rr(n));hold on
%         yyaxis left
%         plot(tBinning,FR2{n},'LineWidth',2,'LineStyle','-');
%         yyaxis right
%         plot(tSingleSti2,SingleSti2,'LineWidth',2,'LineStyle','-')
%     end
%     %============================================
    
%     % All channel raster plot
%     for i=1:60
%         if isempty(Spikes{i})==1
%             Spikes{i}=0;
%         end
%     end
%     figure(z*3);
%     LineFormat.Color = [0 0 0];
%     subplot(5,1,[1:4])
%     plotSpikeRaster(Spikes,'PlotType','vertline','LineFormat',LineFormat)
%     subplot(5,1,5)
%     plot(t,sti)
%     samexaxis('abc','xmt','on','ytac','join','yld',1);
    
    for ii=[1,2]*z
        figure(ii)
        for jj=1:length(On_channel)
            subplot(8,8,rr(On_channel(jj)));hold on
            box on
            ax=gca;
            ax.XColor=[0.8 0 0]; % red
            ax.LineWidth = 1.5;
        end
        for jj=1:length(Off_channel)
            subplot(8,8,rr(Off_channel(jj)));hold on
            box on
            ax=gca;
            ax.XColor=[0 0.6 0]; % green
            ax.LineWidth = 1.5;
        end
        for jj=1:length(OnOff_channel)
            subplot(8,8,rr(OnOff_channel(jj)));hold on
            box on
            ax=gca;
            ax.XColor=[0.6 0 0.6]; % magenta
            ax.LineWidth = 1.5;
        end
    end
end