% Analysis STA experiment, by Rona
clear all
close all
datapath='\\192.168.0.102\Public\Retina\Chou\Exp\20220901\SplitData\';
cd(datapath);
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
SamplingRate = 20000;
cc = hsv(3);
rr =    [9,17,25,33,41,49,...
      2,10,18,26,34,42,50,58,...
      3,11,19,27,35,43,51,59,...
      4,12,20,28,36,44,52,60,...
      5,13,21,29,37,45,53,61,...
      6,14,22,30,38,46,54,62,...
      7,15,23,31,39,47,55,63,...
        16,24,32,40,48,56];
roi = [1:60]; 
color1={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]};
color2={[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0 0]};
mkdir MIandSTA

for z = [16 17 21 22]
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA roi fignames date datapath Date
    file = all_file(z).name;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    bin = 10;  BinningInterval = bin*10^-3;  %ms

%% diode as TriggerData
%     load(['E:\google_rona\20170929\diode\diode_',filename]);
%     [~,locs_a2]=findpeaks(diff(diff(a2)),'MINPEAKHEIGHT',5*std(diff(diff(a2))));
%     TimeStamps_a2 = (locs_a2)/SamplingRate;
%     TriggerData = eyf(TimeStamps_a2(1)*SamplingRate:TimeStamps_a2(end)*SamplingRate);% figure;plot(isi);
%% a_data as TriggerData  
    [b,a] = butter(2,50/20000,'low'); % set butter filter
    a_data2 = filter(b,a,a_data(1,:));
    TriggerData = a_data2(TimeStamps(1)*SamplingRate:TimeStamps(length(TimeStamps))*SamplingRate);% figure;plot(isi);
    inten = downsample(TriggerData,SamplingRate*BinningInterval);
    
    % transform stimulus from volt to intensity
    inten=(inten-32768).*125*10^(-6); % transfrom from raw data of MCRack to output volt
    
    % load calibration data
%     load('\\nasd32ad0\Experiment\Retina\Chou\stimulus saving\19-Apr-2020\calibration\calibration_PAC_19-Apr-2020.mat')
%     inten=inten-offset;
%     Ip=inten/10.421/10^6;
%     r=0.37;
%     P=Ip/r;
%     A=13*10^-6;
%     inten=P/A*1000; % unit: mW/m^2
%% substitute by sorted spikes
%     ss=[29,30,28,27,22,21,14,20,...
%         13,6,12,5,19,11,4,3,10,...
%         18,2,9,1,8,17,7,16,15,26,...
%         25,23,24,32,31,33,34,39,...
%         40,47,41,48,55,49,56,42,...
%         50,57,58,51,43,59,52,60,...
%         53,44,54,45,46,35,36,38,37];
%     Spikes = cell(1,60);
%     
%     load([path_sort,filename(1:end-4),'_sort.mat']);
%     temp_spikes={};
%     for h=1:60
%         if h<11
%             temp_spikes{ss(h)} = eval(['adc00',int2str(h-1)]);
%         else
%             temp_spikes{ss(h)} = eval(['adc0',int2str(h-1)]);
%         end
%     end
%     for i=1:60
%         if isempty(temp_spikes{i})==1
%             continue
%         end
%         for j=1:length(temp_spikes{i}(:,1))
%             if temp_spikes{i}(j,3)==1
%                 Spikes{i}=[Spikes{i} temp_spikes{i}(j,1)];
%             end
%         end
%     end

%================= for sorted data in .xlsx files =========================
%     ss = [29,18,39,53,4,31,59,5,23,38,13,15,...
%         42,37,21,17,55,35,28,9,47,54,10,34,...
%         60,11,32,43,12,25,57,20,16,56,38,22,...
%         8,48,46,30,2,40,44,3,33,52,19,24,...
%         51,6,26,50,14,7,49,36,27,1,41,45];
%     Spikes = {};
%     xls = xlsread([filename(1:end-4),'.xls']);
%     for j = 1:max(xls(:,1))
%             temp = 0;
%         for i = 1:length(xls)
%             if xls(i,1) == j
%                 temp = temp+1;
%                 Spikes{ss(j)}(1,temp) = xls(i,2);
%             end
%         end
%     end
%     BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
%     BinningSpike = zeros(60,length(BinningTime));
%     for i = 1:60
%         [n,xout] = hist(Spikes{i},BinningTime) ;
%         BinningSpike(i,:) = n ;
%     end
%==========================================================================
   
%% spike process
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end
     BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;    
     temprr=0;

%% STA'Name',fignames{z});
    figure(1);
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])

    for nn = 1:length(roi)
        spike = BinningSpike(roi(nn),:);

        window = 1;  %STA window
        window2 = 1;
        sts = [];
        temp = 0;
        spike(1:window/BinningInterval) = 0;
        spike(length(spike)-window2/BinningInterval-10:end) = 0;
        inten2=inten; %-mean(inten)
        for in = 1:length(spike)
           if spike(in)~=0
              temp = temp+1;
              sts(temp,:) = spike(in)*inten2(in-round(window/BinningInterval):in+round(window2/BinningInterval));
           end
        end
        STA = sum(sts)/sum(spike);
        STA = STA/max(abs(STA));
        STA = normalize(STA);
        
        t = [-window*1000:bin:window2*1000];
        subplot(8,8,rr(nn));
        plot(t,STA,'LineWidth',1);%,'color',cc(z,:)
        title(num2str(nn))
        xlim([-window*1000 window2*1000])
%         ylim([-1 1])
        hold on
        STAAAAA{nn}=STA;
        TimeShift{nn}=t;
    end
    
    % ======= plot single channel =============
%     channel=41;
%     figure(2);hold on;box on
%     plot(t,STAAAAA{channel},'LineWidth',2)
%     xlim([-window*1000 window2*1000])
%     ylim([-1 1])
%     xlabel('time (ms)')
%     ylabel('STA')
    
% %======= derivative of  STA ===============
%     for ii=1:60
%         subplot(8,8,rr(ii));hold on
%         yyaxis right
%         if ~isnan(STAAAAA{ii})
%             dSTA=diff(STAAAAA{ii});
%         else
%             dSTA=NaN;
%         end
%         plot(t(1:end-1),dSTA,'LineWidth',1.5)
%     end


% save([datapath,'\MIandSTA\',filename(1:end-4),'_STA.mat'],'STAAAAA','TimeShift')
end



%% Classify ON and OFF cells from cSTA
% onoffindex=[];
% for i =1:60
%     try
%         onoffindex(i)=sum(STAAAAA{1,i})/length(STAAAAA{1,i});
%     catch
%         continue
%     end
% end

%%
% %% change figure color to classify different cells
% % === P ===
% path_class='F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\HMM\';
% file_class=[Date,'_HMM_G30_sort_unit1_MI_also_other_parameters.mat'];
% savepath='F:\我的雲端硬碟\Retina exp\exp data\整理\figures\';
% load([path_class,file_class]);
% figure(1)
% for n=1:length(P_channel)
%     subplot(8,8,rr(P_channel(n)))
%     set(gca,'Color',[0.8 1 0.8])
% end
% % === N ===
% for n=1:length(N_channel)
%     subplot(8,8,rr(N_channel(n)))
%     set(gca,'Color',[0.8 0.8 1])
% end
% 
% % saveas(gcf,[savepath,Date,'_STA_NPclassified.fig']);
%    