% find different cell types form Spike Trigger
% Analysis STA experiment, by Rona
clear all
close all
Date='20200419';
cd(['F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\',Date,'\']);
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

for z = [76]
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA roi fignames date path_sort Date
    file = all_file(z).name ;
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

%% spike process
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end
     BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;
     temprr=0;


    for nn = 39%1:length(roi)
        %% STA
        spike = BinningSpike(roi(nn),:);

        window = 1;  %STA window
        window2 = 1;
        sts = [];
        temp = 0;
        spike(1:window/BinningInterval) = 0;
        spike(length(spike)-window2/BinningInterval-10:end) = 0;
        inten2=inten-mean(inten);
        for in = 1:length(spike)
           if spike(in)~=0
              temp = temp+1;
              sts(temp,:) = spike(in)*inten2(in-round(window/BinningInterval):in+round(window2/BinningInterval));
           end
        end
        STA = sum(sts)/sum(spike);
        STAstd=std(sts);
        [MSTA,iM]=max(STA);
        [mSTA,im]=min(STA);
        
        figure(1);hold on
        t = [-window*1000:bin:window2*1000];
        plot(t,STA,'-',t,STAstd/2+STA,'--',t,-STAstd/2+STA,'--','color',[0 0.4470 0.7410],'LineWidth',1.5);%,'color',cc(z,:)
        plot(t(iM),MSTA,'o')
        plot(t(im),mSTA,'o')
        xlim([-window*1000 window2*1000])
%         ylim([-1 1])
        
        for i=1:size(sts,1)
            st1(i)=sts(i,im); % min peak value of ST
        end
        [stoff,noff]=find(st1<0);
        
        for i=1:size(sts,1)
            st2(i)=sts(i,iM); % max peak value of ST
        end
%         ntot=1:size(sts,1);
%         nttt=ntot;
%         for i=ntot
%             nttt(noff)=0;
%         end
%         noff_other=find(nttt);
        [ston,non]=find(st2>0); %(noff_other)
%         non=noff_other(non);
        
        [C,ia,ib] = intersect(noff,non);
        n_biphasic=noff(ia);
        [non1,ii]=setdiff(non,n_biphasic);
        [noff1,jj]=setdiff(noff,n_biphasic);
        
        
        figure(2)
        try
            plot(t,sts(noff1,:))
            title('off ST')
        catch
        end
        figure(3)
        plot(t,sts(non1,:))
        title('on ST')
        figure(4)
        plot(t,sts(n_biphasic,:))
        title('on off ST')
        
        
        
        STA2=sum(sts(noff1,:))/length(noff);
        STAstd2=std(sts(noff1,:));
        STA3=sum(sts(non1,:))/length(non);
        STAstd3=std(sts(non1,:));
        STA4=sum(sts(n_biphasic,:))/length(n_biphasic);
        STAstd4=std(sts(n_biphasic,:));
        
        figure;
        subplot(1,3,1);plot(t,STA2,'-',t,STAstd2/2+STA2,'--',t,-STAstd2/2+STA2,'--','LineWidth',2,'color',[0 0.4470 0.7410])    
        title(['Off type STA  ',num2str(round(length(noff1)/size(sts,1)*100)),'%'],'fontsize',20)
        xlabel('time shift (ms)','fontsize',16)
        ylabel('intensity (not absolute value)','fontsize',16)
%         ylim([-150 150])
        subplot(1,3,2);plot(t,STA3,'-',t,STAstd3/2+STA3,'--',t,-STAstd3/2+STA3,'--','LineWidth',2,'color',[0.8500 0.3250 0.0980])
        title(['On type STA  ',num2str(round(length(non1)/size(sts,1)*100)),'%'],'fontsize',20)
        xlabel('time shift (ms)','fontsize',16)
%         ylim([-150 150])
        subplot(1,3,3);plot(t,STA4,'-',t,STAstd4/2+STA4,'--',t,-STAstd4/2+STA4,'--','LineWidth',2,'color',[0.9290 0.6940 0.1250])
        title(['On-Off type STA  ',num2str(round(length(n_biphasic)/size(sts,1)*100)),'%'],'fontsize',20)
        xlabel('time shift (ms)','fontsize',16)
%         ylim([-150 150])
        
        
        disp(['biphasic Spike Trigger: ',num2str(length(n_biphasic)),'  ',num2str(length(n_biphasic)/size(sts,1)*100),'%'])
        disp(['On Spike Trigger: ',num2str(length(non1)),'  ',num2str(length(non1)/size(sts,1)*100),'%'])
        disp(['Off Spike Trigger: ',num2str(length(noff1)),'  ',num2str(length(noff1)/size(sts,1)*100),'%'])
        disp(['total Spike Triggers: ',num2str(size(sts,1))])
        
        %% MI
        
        
    end
end