% Analysis STA experiment, by Rona
clear all
close all
datapath='\\192.168.0.102\Public\Retina\Chou\Exp\20220916\SplitData';
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
file_numbers=[31 33 32 34 35]; % file selections 16 18 17 19 20;31 33 32 34 35;21 23 22 29 30

colors_default={[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560], ...
    [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};
color1={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]};
color2={[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0 0]};
ls_set={'-','--',':','-.'};

mkdir MIandSTA

% for i_set = 1:size(file_numbers,1)
for z = 1:length(file_numbers)
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr...
        information rr STA roi fignames date datapath Date colors_default file_numbers ls_set i_set BinningSpike_seq...
        nonlinear_in nonlinear_out kt_set
    file = all_file(file_numbers(z)).name;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    bin = 25;  BinningInterval = bin*10^-3;  %ms

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
%   load calibration data
    load('\\192.168.0.102\Public\Retina\Chou\Exp\20220916\16-Sep-2022\calibration\calibration_PAC_16-Sep-2022.mat')
    inten=inten-offset;
    Ip=inten/10.421/10^6;
    r=0.37;
    P=Ip/r;
    A=13*10^-6;
    inten=P/A*1000; % unit: mW/m^2
   
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
    
        window = 2;  %STA window
        window2 = 0;
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
        
        STA = sum(sts,1)/sum(spike);
    
%         STA = STA/max(abs(STA));
%         STA = normalize(STA);
        
        t = [-window*1000:bin:window2*1000];
        sub1=subplot(8,8,rr(nn));
        try
            plot(t,STA,'LineWidth',1);% ,'color',colors_default{z} ,'linestyle',ls_set{i_set}
        catch
        end
        title(num2str(nn))
        xlim([-window*1000 window2*1000])
    %         ylim([-1 1])
        hold on
        STAAAAA{nn}=STA;
        TimeShift{nn}=t;
    end
    
    %% calculate convoluted signal and obtain nonlinear filter
    channel=22;
    kt_set{z}=STAAAAA{channel}-mean(STAAAAA{channel});
%     kt_set{z}=normalize(STAAAAA{channel});
    gt=conv(inten,kt_set{z}(end:-1:1),'valid'); % (length(kt)-1) elements are dropped
    figure(2);hold on
    nonlinear_in(z,:)=gt(1:end);
    nonlinear_out(z,:)=BinningSpike(channel,length(kt_set{z}):end);
    plot(nonlinear_in(z,:),nonlinear_out(z,:)+(z-1)*0.1,'o')
    figure(3)
    subplot(2,1,1);hold on;plot(kt_set{z})
    subplot(2,1,2);hold on;plot(gt)
    
end

%% find mean curve of nonlinearity
std_gt_value=zeros(length(file_numbers),max(nonlinear_out,[],'all'));
mean_gt_value=zeros(length(file_numbers),max(nonlinear_out,[],'all'));
for nfile=1:length(file_numbers)
    for n_spike=1:max(nonlinear_out(nfile,:))
        inds=find(nonlinear_out(nfile,:)==n_spike);
        if isempty(inds)
        else
            std_gt_value(nfile,n_spike)=std(nonlinear_in(nfile,inds));
            mean_gt_value(nfile,n_spike)=mean(nonlinear_in(nfile,inds));
        end
   
    end
    figure(4);hold on
    errorbar(mean_gt_value(nfile,:),1:size(mean_gt_value,2)+nfile*0.1,std_gt_value(nfile,:),'horizontal','o','linewidth',2)
end
%% standardize the nonlinearity
i=0;
n_zero=0;
while n_zero<1
    i=i+1;
    ind_0=find(mean_gt_value(:,i)==0);
    n_zero=length(ind_0);
    useful_states=i-1;
end
figure(5);subplot(1,2,1);hold on
fraction_NL=zeros(size(mean_gt_value,1)-1,useful_states);
mean_frac_NL=ones(size(mean_gt_value,1),1);
for j=1:size(mean_gt_value,1)-1
    fraction_NL(j,:)=mean_gt_value(end,1:useful_states)./mean_gt_value(j,1:useful_states);
    mean_frac_NL(j)=mean(fraction_NL(j,:));
    plot(mean_gt_value(j,:)*mean_frac_NL(j),1:size(mean_gt_value,2),'o','linewidth',2)
end
plot(mean_gt_value(end,:),1:size(mean_gt_value,2),'o','linewidth',2)
subplot(1,2,2);hold on
for i=1:size(kt_set,2)
    plot(kt_set{i}*mean_frac_NL(i))
end

