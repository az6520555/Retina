% NP and P cell same spike analysis
close all;
clear all;

channel_number = [58 39]; %2 element only for now. 1st for N, 2nd for P.
np='NP';
tolerance = 0.01; %s
name='OU_tau=600ms_19-Apr-2020_0_sort_unit1'; %_cutoff=10

load(['F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200419\',name,'.mat']) % load spike time file
load(['F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200419\MIandSTA\',name,'_MI.mat']) % load original mutual information

% t=0:1/rate:(size(a_data,2)-1)/rate;
% for i=1:60
%     Spikes{i}=Spikes{i}(Spikes{i}>TimeStamps(1) & Spikes{i}<TimeStamps(2))-TimeStamps(1);
% end
analyze_spikes = Spikes;
bin=10;  BinningInterval = bin*10^-3; 
SamplingRate=20000;

%% plot original MI
figure(1);
for k = channel_number
    plot(TimeShift,MI{k},'LineWidth',1.5,'LineStyle','-');hold on;
    %plot(time,smooth(Mutual_shuffle_infos{j}),'LineWidth',1.5,'LineStyle','-');hold on;
    xlim([-1000 1000])
%     ylim([0 inf+0.5])
    xlabel('\deltat (ms)');ylabel('MI(\gamma(t),I(t-\deltat)) (bits/s)');
    set(gca,'fontsize',12); hold on
    title('MI(\gamma(t),I(t-\deltat)) of N and P type cell');
end 
legend('N type cell','P type cell')

%%  compare the raster plot of different cells
% for i=1:60
%     if isempty(Spikes{i})==1
%         Spikes{i}=0;
%     end
% end
% figure(87);
% LineFormat.Color = [0.3 0.3 0.3];
% plotSpikeRaster({Spikes{26} Spikes{40}},'PlotType','vertline','LineFormat',LineFormat)
% hold on

%% remove sharing spike and plot the cancelled, sharing, and original spike
for j = 1:length(channel_number)
    sub_Spikes{j} = analyze_spikes{channel_number(j)};
    if isempty(sub_Spikes{j})==1
        sub_Spikes{j}=0;
    end;
    
end
null_index_1 = [];
null_index_2 = [];
%There may be a better way to search
for m = 1:length(sub_Spikes{1})
    for n = 1:length(sub_Spikes{2})
        if abs(sub_Spikes{1}(m)-sub_Spikes{2}(n)) < tolerance
            null_index_1 = [null_index_1 m];
            null_index_2 = [null_index_2 n];
        end
    end
end
null_index_1 = unique(null_index_1);
null_index_2 = unique(null_index_2);
sub_Spikes{3} = sub_Spikes{1}(null_index_1);
sub_Spikes{4} = sub_Spikes{2}(null_index_2);
sub_Spikes{5} = sub_Spikes{1};
sub_Spikes{6} = sub_Spikes{2};
sub_Spikes{1}(null_index_1) = [];
sub_Spikes{2}(null_index_2) = [];
figure(2);subplot(5,1,[1,4])
LineFormat.Color = [0.1 0.1 0.7];
plotSpikeRaster({sub_Spikes{5} 0 sub_Spikes{3} 0 sub_Spikes{1} 0},'PlotType','vertline','LineFormat',LineFormat);hold on
LineFormat.Color = [0.7 0.1 0.1];
plotSpikeRaster({0 sub_Spikes{6} 0 sub_Spikes{4} 0 sub_Spikes{2}},'PlotType','vertline','LineFormat',LineFormat);

hold on;
oN = length(analyze_spikes{channel_number(1)}); % total number of spikes in cell 1 
oP = length(analyze_spikes{channel_number(2)}); % total number of spikes in cell 2
analyze_spikes{channel_number(1)}(null_index_1) = [];
analyze_spikes{channel_number(2)}(null_index_2) = [];
cN = length(analyze_spikes{channel_number(1)}); % total number of spikes which common spikes are cancelled in cell 1
cP = length(analyze_spikes{channel_number(2)}); % total number of spikes which common spikes are cancelled in cell 2
disp(['independent spikes ratio of cell 1: ',num2str(cN/oN)])
disp(['independent spikes ratio of cell 2: ',num2str(cP/oP)])
% ============
% subspikes{1} independent spikes of cell 1
% subspikes{2} independent spikes of cell 2
% subspikes{3} matched spikes of cell 1 
% subspikes{4} matched spikes of cell 2
% subspikes{5} original spikes of cell 1
% subspikes{6} original spikes of cell 2

%% stimulus
[b,a] = butter(2,50/20000,'low'); % set butter filter
a_data2 = filter(b,a,a_data(1,:));
isi = a_data2(TimeStamps(1)*20000:TimeStamps(length(TimeStamps))*20000); % filtered stimulus data

T_isi=1/SamplingRate:1/SamplingRate:length(isi)/SamplingRate;
T_isi=T_isi+TimeStamps(1);
figure(2);subplot(5,1,5)
plot(T_isi,isi);
samexaxis('abc','xmt','on','ytac','join','yld',1)

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

%% Spikes binning process

BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
BinningSpike = zeros(length(sub_Spikes),length(BinningTime));
for i = 1:length(sub_Spikes)
    [n,xout] = hist(sub_Spikes{i},BinningTime);
    BinningSpike(i,:) = n;
end
BinningSpike(:,1) = 0; BinningSpike(:,end) = 0;

%% calculate and plot  MI of cancelled spike
ana_train=[5 6 4 1 2];

linestyle={'-','-','--',':',':'};

backward=ceil(1000/bin); forward=ceil(1000/bin);
figure('Name','MI of different data');hold on;box on
for nn= 1:length(ana_train) %1:length(sub_Spikes)
    Neurons = BinningSpike(ana_train(nn),:); 
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
%     corr=[corrp corrf];
%     fr(z)=mean(Neurons)/BinningInterval; %firing rate(Hz)
[pks(nn),plocs(nn)]=max(information);
t=[-backward*bin:bin:forward*bin];  

plot(t,information,'LineWidth',2,'LineStyle',linestyle{nn});%,'color',cc(z,:)

end
% title('MI(\gamma(t),I(t-\deltat)) of original P cell and matched P cell')
xlabel('\deltat (ms)');ylabel('bits/s'); %MI(\gamma(t),I(t-\deltat)) 
% legend(['P'],['P^{\prime}'])
% legend(['Original N'],['Original P'], ['Cancelled N'],['Cancelled P'], ['sharing N'],['sharing P']);
legend('N','P','P^{\prime}','N^{\prime}^{\prime}','P^{\prime}^{\prime}')

%% STA

BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
BinningSpike = zeros(2,length(BinningTime));
for i=[1 2]
    sspi_ind=[4 2];
    [n,xout] = hist(sub_Spikes{sspi_ind(i)},BinningTime) ;
    BinningSpike(i,:) = n ;
end
BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;    
temprr=0;


[b,a] = butter(2,50/20000,'low'); % set butter filter
a_data2 = filter(b,a,a_data(1,:));
TriggerData = a_data2(TimeStamps(1)*SamplingRate:TimeStamps(length(TimeStamps))*SamplingRate);% figure;plot(isi);
inten = downsample(TriggerData,SamplingRate*BinningInterval);

xxx={'-','--'};
for nn = 1:size(BinningSpike,1)
    spike = BinningSpike(nn,:);

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

%         figure(roi(nn));hold on;
    figure(4);hold on
    t = [-window*1000:bin:window2*1000];    
    plot(t,STA,'LineWidth',1.5,'linestyle',xxx{nn},'color','r');%,'color',cc(z,:)
    xlim([-window*1000 window2*1000])
    hold on
    ylabel('STA')
    legend('P^{\prime}','P^{\prime}^{\prime}')

%     shift_bin=10;
%     t_STA=-1000:shift_bin:1000;t_STA=t_STA(1:end-1)+shift_bin/2;
%     if ~isnan(STA)
%         dSTA=diff(STA);
%     else
%         dSTA=NaN;
%     end
%     dSTA=dSTA/max(abs(dSTA));

%     figure(nn*10);hold on
%     t = [-window*1000:bin:window2*1000];    
%     yyaxis right
%     plot(t_STA,dSTA,'LineWidth',2);%,'color',cc(z,:)
%     xlim([-window*1000 window2*1000])
%     hold on
%     ylabel('dSTA/dt')
end
xlabel('time before and after the spike')
% title('STA of original P cell and matched P cell')
% legend(['STA of original P'],['STA of matched P'])

%% calculate mutual information between intensity derivative and firing rate

  % state of changing rate
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


backward=ceil(3000/bin); forward=ceil(3000/bin);
for nn= 1:length(ana_train)
    Neurons = BinningSpike(nn,:); 
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
    %     corr=[corrp corrf];
    %     fr(z)=mean(Neurons)/BinningInterval; %firing rate(Hz)
    [pks(nn),plocs(nn)]=max(information);
    t=[-backward*bin:bin:forward*bin];
%     figure(nn*10);hold on
%     yyaxis left
%     plot(t,information,'LineWidth',2,'LineStyle','-');%,'color',cc(z,:)
%     xlabel('\deltat (ms)');ylabel('MI(\gamma(t),dI(t-\deltat)/dt) (bits/s)');
end




