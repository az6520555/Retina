% plot MI and STA at the same time
clear all
close all

Date='20190124';
path=['F:\我的雲端硬碟\Retina exp\exp data\Sorted_final_data\',Date];
cd(path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file);
SamplingRate=20000;
file_seq=[7];
corr_time=[40 160 300 540];
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
color1={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]};
color2={[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0 0]};
neg=1.5;
pos=1.5;

mkdir MIdata

for z = 1:length(file_seq)
    file = all_file(file_seq(z)).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    bin=10;  BinningInterval = bin*10^-3; 
    backward=ceil(neg*1000/bin); forward=ceil(pos*1000/bin);
    
    
    [b,a] = butter(2,50/20000,'low'); % set butter filter
    a_data2 = filter(b,a,a_data(3,:));
    isi = a_data2(TimeStamps(1)*20000:TimeStamps(length(TimeStamps))*20000);% figure;plot(isi);
    
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        if length(Spikes{i})<100
            Spikes{i}=[];
        end
        [n,xout] = hist(Spikes{i},BinningTime);
        BinningSpike(i,:) = n;
    end
    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;% figure;plot(BinningTime,sum(BinningSpike),BinningTime,10*Stimuli,'o')

    
       %% state of light intensity %%% 
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
    
    MI = cell(1,60); % create an array to save the MI data
    xcor = cell(1,60); % create an array to save the cross correlation data
    infor=[];co=[];
    for nn= roi
        n=roi(nn);
        Neurons = BinningSpike(n,:); 
        
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
    
    % =================multifigure==========================
%     figure(1);subplot(8,8,rr(n));hold on
%     yyaxis left
%     plot(t,information,'LineWidth',2,'LineStyle','-','color',color1{z});%,'color',cc(z,:)
    % ==================================================
    
%     yyaxis right
%     plot(t,corr,'LineWidth',2,'LineStyle','-')
%     refline([0 0])
%     xlim([-1000 1000])
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
    MI{nn} = information;
    xcor{nn} = corr;
    
    end
    
    %% STA'Name',fignames{z});
    for nn = roi
        spike = BinningSpike(roi(nn),:);

        window = neg;  %STA window
        window2 = pos;
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
        STA = sum(sts)/sum(spike);%;
        STA = STA/max(abs(STA));
        
%=======================multifigure============================
%         subplot(8,8,rr(nn));
%         t = [-window*1000:bin:window2*1000];
%         yyaxis right
%         plot(t,STA,'LineWidth',2,'LineStyle','-','color',color2{z});%,'color',cc(z,:)
%         refline([0 0])
%         xlim([-window*1000 window2*1000])
%         ylim([-1 1])
%         hold on
%===================================================
        
        STAAAAA{nn}=STA;
    end
    
    %======================== single plot=============================
    channel=39;
    figure(10);
%     subplot(2,2,z);
    hold on
    yyaxis left
    plot(t,MI{channel},'LineWidth',2);%,'LineStyle','-','color',color1{z}
    ylabel('mutual information (bits/s)','FontSize',16)
    yyaxis right
    plot(t,STAAAAA{channel},'LineWidth',2)%,'LineStyle',':','color',color1{z}
    ylabel('STA','FontSize',16)
    xlabel('time shift (ms)','FontSize',16)
    refline([0 0])
    xlim([-neg*1000 pos*1000])
    title(['correlation time = ',num2str(corr_time(z)),'ms'],'FontSize',20)
    %=====================================================================

    
end

%% single plot
% channel=26;
% figure
% yyaxis left
% plot(t,MI{channel},'LineWidth',2,'LineStyle','-');%,'color',cc(z,:)
% ylabel('mutual information (bits/s)')
% yyaxis right
% plot(t,STA,'LineWidth',2,'LineStyle','-')
% ylabel('STA')
% xlabel('time shift (s)')
% refline([0 0])
% xlim([-1000 1000])

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

% saveas(gcf,[savepath,Date,'_STA_NPclassified.fig']);
   