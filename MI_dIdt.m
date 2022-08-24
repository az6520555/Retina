% Calculated MI for continuous changing intensity stimulation , by Rona
% MI between intensity changing rate and firing rate

clear all
close all



Date='20190408';
path=['F:\我的雲端硬碟\Retina exp\exp data\Sorted_final_data\',Date];
MI_inten_file='F:\我的雲端硬碟\Retina exp\exp data\Sorted_final_data\20190131\MIdata\20190131_HMM_G5_sort_unit1_MI.mat'; %
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
roi = [1:60];
corr_time=[40 160 300 540];

mkdir MIdata
file_seq=[7 1 10 4];
channel=51;

for z = 1:length(file_seq)
    file = all_file(file_seq(z)).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    bin=10;  BinningInterval = bin*10^-3; 
    backward=ceil(1000/bin); forward=ceil(1000/bin);
%%  TimeStamps  %%%%%%%%%%%%%%%%%%%
%     if size(a_data,1)==1              %only find one analog channel, possibly cause by the setting in MC_rack
%         a_data2 = a_data(1,:);
%     else
%         a_data2 = a_data(2,:);   
%     end
%     [~,locs]=findpeaks(diff(a_data2),'MINPEAKHEIGHT',5*std(diff(a_data2)));
%     analog_loc = (locs)/SamplingRate;
%     TimeStamps = analog_loc;
%    
% TimeStamps =[TimeStamps_H(1):TimeStamps_N(end)]
    if length(TimeStamps)==1
        TimeStamps(2)=TimeStamps(1)+200;
    end

   %% diode as isi2
%    load(['E:\google_rona\20170502\diode\diode_',filename]);
%      [~,locs_a2]=findpeaks(diff(diff(a2)),'MINPEAKHEIGHT',5*std(diff(diff(a2))));
%      TimeStamps_a2 = (locs_a2)/SamplingRate;
% 
% %     [b,a]=butter(2,60/20000,'low');
% %     a_data2=filter(b,a,callumin)';
% %     a_data2=eyf;
%     isi = callumin_filter(TimeStamps_a2(1)*20000:TimeStamps_a2(end)*20000);
% %     figure(z);autocorr(isi,100000);
    %% a_data as isi  %   a_data2=(a_data-32768)*0.1042;
    [b,a] = butter(2,30/20000,'low'); % set butter filter
    a_data2 = filter(b,a,a_data(3,:));
    isi = a_data2(TimeStamps(1)*20000:TimeStamps(length(TimeStamps))*20000);% figure;plot(isi);
 %% Spike process
 % if the data is sorted, you can load the sorted excel file here
%    ss = [29,18,39,53,4,31,59,5,23,38,13,15,...
%     42,37,21,17,55,35,28,9,47,54,10,34,...
%     60,11,32,43,12,25,57,20,16,56,38,22,...
%     8,48,46,30,2,40,44,3,33,52,19,24,...
%     51,6,26,50,14,7,49,36,27,1,41,45];
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
 
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        if length(Spikes{i})<100
            Spikes{i}=[];
        end
        [n,xout] = hist(Spikes{i},BinningTime);
        BinningSpike(i,:) = n;
    end
    % [n,out] = hist(TimeStamps,BinningTime);
    % Stimuli = n;
    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;% figure;plot(BinningTime,sum(BinningSpike),BinningTime,10*Stimuli,'o')
%     figure;imagesc(BinningTime,[1:60],BinningSpike)
    
    %% state of light intensity changing %%% 
    isi2=[];
    states=8;
    X=isi;
    nX = sort(X);
    abin = length(nX)/states;
    intervals = [nX(1:abin:end) inf]; 
    temp=0;
    for jj = 1:BinningInterval*SamplingRate:length(X)
        temp=temp+1;
        isi2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
        inten(temp)= X(jj);
    end
%     figure(z*10);autocorr(inten,100)
% %     figure;hist(isi2,[1:states]);
% %     figure(50);plot(isi2,'color',cc(z,:));hold on
%% timewindow of stim %%%
%     isi = TimeStamps(2:1:end) - TimeStamps(1:1:end-1); isi = isi*1000;  %for pseudo-period
%     X = zeros(size(BinningTime));
%     temp = 1;
%     for ii = 1:length(X)
%         if BinningTime(ii)<=TimeStamps(temp)
%             X(ii) = isi(temp);%find(isi<TimeStamps(temp),1);
%         else
%             temp = temp + 1;
%             if temp>length(isi); break; end
%             X(ii) = isi(temp);
%         end
%     end
% 
%     m = mean(isi);  %% tt = abs(isi-m);  dev = max(tt(find(tt<300)));
%     isi2=[];
%     states=25;
%     nX = sort(X);
%     abin = length(nX)/states;
%     intervals = [nX(1:abin:end) inf]; 
%     for jj = 1:length(X)
%         [a,b] = find(X(jj)<intervals,1); % stimulus for every 50ms
%         isi2(jj) = b-1;
%     end
  %% state of changing rate %%% 
    isi3=[];isi4=[];
    temp=0;
    for j=1:BinningInterval*SamplingRate:length(isi)
        temp=temp+1;
        isi3(temp) = isi(j);
    end
    
    for i=1:length(isi3)-1
        isi4(i)=isi3(i+1)-isi3(i);
    end
    isi4 = smooth(isi4,10)';
    isi2=[];
    states=8;
    X=isi4;
    nX = sort(X);
    abin = length(nX)/states;
    intervals = [nX(1:abin:end) inf]; 
    temp=0;
    for jj = 1:length(X)
        temp=temp+1;
        isi2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
    end

%     figure;hist(isi2,[1:states]);

%% Mutual Information
MI = cell(1,60); % create an array to save the MI data
xcor = cell(1,60); % create an array to save the cross correlation data
infor=[];co=[];
    for nn= roi
        n=nn;
        Neurons = BinningSpike(n,:); 
%         Neurons = sum(BinningSpike(:,:));
         %% shuffle
%         r=randperm(length(Neurons));
%         for j=1:length(r)            
%             sNeurons(j)=Neurons(r(j));
%         end
%         Neurons=sNeurons;

    dat=[];informationp=[];temp=backward+2;
    for i=1:backward+1 %past(t<0)
        temp=temp-1;
        if length(Neurons)==length(isi2)
            x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
        else
            x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1)-1)';
        end    
        y=isi2(forward+1:length(isi2)-backward)';
        dat=[x,y];
        [N,C]=hist3(dat,[max(Neurons)+1,8]); %20:dividing firing rate  6:# of stim
        px=sum(N,1)/sum(sum(N)); % x:stim
        py=sum(N,2)/sum(sum(N)); % y:word
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
              temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
            end
        end
       
        informationp(temp)=nansum(temp2(:));
        c=corrcoef(x,y);
        corrp(temp)=c(2,1);
%         informationp(temp) = ImExtrapolation_function(dat,8);
    end  
    
    dat=[];informationf=[];temp=0;sdat=[];
    for i=1:forward
        temp=temp+1;
        if length(Neurons)==length(isi2)
            x =Neurons(forward+1-i:length(Neurons)-backward-i)';
        else
            x = Neurons(forward+1-i:length(Neurons)-backward-i-1)';
        end     
        y = isi2(forward+1:length(isi2)-backward)';
        dat=[x,y];
        [N,C]=hist3(dat,[max(Neurons)+1,8]); %20:dividing firing rate  6:# of stim
        px=sum(N,1)/sum(sum(N)); % x:stim
        py=sum(N,2)/sum(sum(N)); % y:word
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
                temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
            end
        end
        
        informationf(temp)=nansum(temp2(:)); 
%      
        c=corrcoef(x,y);
        corrf(temp)=c(2,1);
%             informationf(temp) = ImExtrapolation_function(dat,8);
    end

    information=[informationp informationf]/BinningInterval;
    corr=[corrp corrf];
%     fr(z)=mean(Neurons)/BinningInterval; %firing rate(Hz)
    [pks(z,nn),plocs(z,nn)]=max(information);
    t=[-backward*bin:bin:forward*bin];  
    
    %==========================multifigure===================================
%     figure(1);subplot(8,8,rr(n));hold on;n
% %     yyaxis left
%     plot(t,information,'LineWidth',2,'LineStyle','-');%,'color',cc(z,:)
%     plot([0 0], ylim, '-r')
%     xlim([-1000 1000])
    %======================================================================
    
%     
%     figure(2);hold on
%     plot(t,information,'LineWidth',2,'LineStyle','-');%,'color',cc(z,:)
%     ylabel('mutual information (bits/s)')
%     xlabel('time shift (s)')
%     xlim([-1000 1000])
%     plot([0 0], ylim, '-r')
%     ax = gca;
%     ax.XGrid = 'off';
%     ax.YGrid = 'on';
    
    % ==================
    dMI{nn} = information;
    TimeShift=t;
    save([path,'\MIdata\',filename(1:end-4),'_dMI.mat'],'dMI','TimeShift')
    
    end


        %% STA'Name',fignames{z});
    for nn = roi
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
        STA = STA/max(abs(STA));
        
%=======================multifigure============================
%         subplot(8,8,rr(nn));
%         t = [-window*1000:bin:window2*1000];
%         yyaxis right
%         plot(t,STA,'LineWidth',2,'LineStyle','-');%,'color',cc(z,:)
%         refline([0 0])
%         xlim([-window*1000 window2*1000])
%         ylim([-1 1])
%         hold on
%===================================================
        
        STAAAAA{nn}=STA;
    end
    
    shift_bin=10;
    t_STA=-1000:shift_bin:1000;t_STA=t_STA(1:end-1)+shift_bin/2;
    for ii=roi
        
        if ~isnan(STAAAAA{ii})
            dSTA=diff(STAAAAA{ii});
        else
            dSTA=NaN;
        end
%         dSTA=dSTA/max(abs(dSTA));
        dSTAAAAA{ii}=dSTA;
        %=========================multifigure==============================
%         subplot(8,8,rr(ii));hold on
%         yyaxis right
%         plot(t_STA,dSTA,'LineWidth',2,'LineStyle','-')
%         refline([0 0])
%         ylim([-1 1])
        %==================================================================
    end
    
    %% plot figure independently
    
    figure(2)
    subplot(2,2,z)
    title(['correlation time = ',num2str(corr_time(z))],'FontSize',20)
    yyaxis left
    plot(t,dMI{channel},'LineWidth',2,'LineStyle','-');
    ylabel('MI_\gamma_,_d_I_/_d_t','FontSize',18) % MI[\gamma(t),dI(t-\deltat)/dt]
    yyaxis right 
    plot(t_STA,dSTAAAAA{channel},'LineWidth',2,'LineStyle','-')
    ylabel('$\frac{d}{dt}STA$','Interpreter','latex','FontSize',18)
    xlabel('time shift (ms)','FontSize',18)
    refline([0 0])
%     ylim([-1 1])
    xlim([-1000 1000])
    
    
end

%% plot MI dIdt peak timeshift
% figure(1)
% for i=roi
%     subplot(8,8,rr(i));hold on
%     [dmipk{i},dmilocs{i}]=findpeaks(dMI{i},'MinPeakDistance',10,'MinPeakHeight',max(dMI{i})/2);
% %     if dmipk(1)*dmipk(2)>0
%         
%     yyaxis right
%     if ~isnan(dSTAAAAA{i})
%         plot(t_STA(dmilocs{i}),dSTAAAAA{i}(dmilocs{i}),'ko')
%     else
%         continue
%     end
% end
%     


%% plot MI intensity 
% figure(1)
% for i=roi
%     load(MI_inten_file)
%     subplot(8,8,rr(i));
%     yyaxis left
%     plot(TimeShift,MI{i},'LineStyle','-','LineWidth',2,'color',[0.6350 0.0780 0.1840]);%,'color',cc(z,:)
%     xlim([-1000 1000])
%     hold on
% end

%% change figure color to classify different cells
% === P ===
path_class='F:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\HMM\';
file_class=[Date,'_HMM_G30_sort_unit1_MI_also_other_parameters.mat'];
savepath='F:\我的雲端硬碟\Retina exp\exp data\整理\figures\';
load([path_class,file_class]);
figure(1)
for n=1:length(P_channel)
    subplot(8,8,rr(P_channel(n)))
    set(gca,'Color',[0.8 1 0.8])
end
% === N ===
for n=1:length(N_channel)
    subplot(8,8,rr(N_channel(n)))
    set(gca,'Color',[0.8 0.8 1])
end





