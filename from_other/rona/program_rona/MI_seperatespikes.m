%% Seperate spike into common and uncommon spike between P and N 
clear all
close all
date='20181028';
cd(['E:\Chou\',date]);

all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
SamplingRate=20000;
cc=hsv(n_file);

for z = 5
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA ro date
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    bin=10;  BinningInterval = bin*10^-3; 
%%  TimeStamps  %%%%%%%%%%%%%%%%%%%
    if size(a_data,1)==1              %only find one analog channel, possibly cause by the setting in MC_rack
        a_data2 = a_data(1,:);
    else
        a_data2 = a_data(2,:);   
    end
    [~,locs]=findpeaks(diff(a_data2),'MINPEAKHEIGHT',5*std(diff(a_data2)));
    analog_loc = (locs)/SamplingRate;
    TimeStamps = analog_loc;
%     TimeStamps = [30,300];
    if length(TimeStamps)==1
        TimeStamps(2)=TimeStamps(1)+200;
    end

   %% diode as isi2
%    load(['E:\google_rona\20161115\diode\diode_',filename]);
%      [~,locs_a2]=findpeaks(diff(diff(a2)),'MINPEAKHEIGHT',5*std(diff(diff(a2))));
%      TimeStamps_a2 = (locs_a2)/SamplingRate;
%      
%     [b,a]=butter(2,60/20000,'low');
%     a_data2=filter(b,a,callumin)';
% %     a_data2=eyf;
%     isi = a_data2(TimeStamps_a2(1)*20000:TimeStamps_a2(end)*20000);
    %% a_data as isi 
% %   a_data2=(a_data-32768)*0.1042;
    a_data2=a_data(3,:);
    isi = a_data2(TimeStamps(1)*20000:TimeStamps(end)*20000);% figure;plot(isi);
    m = mean(isi); % mV
    DataTime =TimeStamps(end)-TimeStamps(1); %s
    
    %% Spike process
    % sorting data
%      xls=xlsread([name,'.xls']);
%         ss = [29,18,39,53,4,31,59,5,23,38,13,15,...
%         42,37,21,17,55,35,28,9,47,54,10,34,...
%         60,11,32,43,12,25,57,20,16,56,38,22,...
%         8,48,46,30,2,40,44,3,33,52,19,24,...
%         51,6,26,50,14,7,49,36,27,1,41,45];
%     Spikes=cell(1,60);
%     for j = 1:max(xls(:,1))
%             temp = 0;
%         for i = 1:length(xls)
%             if xls(i,1) == j
%                 temp = temp+1;
%                 Spikes{ss(j)}(1,temp) = xls(i,2);
%             end
%         end
%     end
% =============================================================
%     ss = [29,18,39,53,4,31,59,5,23,38,13,15,...
%         42,37,21,17,55,35,28,9,47,54,10,34,...
%         60,11,32,43,12,25,57,20,16,56,38,22,...
%         8,48,46,30,2,40,44,3,33,52,19,24,...
%         51,6,26,50,14,7,49,36,27,1,41,45];
%     Spikes = cell(1,60);
%     path_sort=['I:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\sorted data\',date,'\'];
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
% =================================================================

    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end
    % [n,out] = hist(TimeStamps,BinningTime);
    % Stimuli = n;
    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;% figure;plot(BinningTime,sum(BinningSpike),BinningTime,10*Stimuli,'o')
%     figure;imagesc(BinningTime,[1:60],BinningSpike)
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
% %     figure;hist(isi2,[1:states]);
% %     figure(50);plot(isi2,'color',cc(z,:));hold on

%% state of changing rate %%% 
%     isi3=[];isi4=[];
%     temp=0;
%     for j=1:BinningInterval*SamplingRate:length(isi)
%         temp=temp+1;
%         isi3(temp) = isi(j);
%     end
%     
%     for i=1:length(isi3)-1
%         isi4(i)=isi3(i+1)-isi3(i);
%     end
%     isi4 = smooth(isi4,5)';
%     isi2=[];
%     states=8;
%     X=isi4;
%     nX = sort(X);
%     abin = length(nX)/states;
%     intervals = [nX(1:abin:end) inf]; 
%     temp=0;
%     for jj = 1:length(X)
%         temp=temp+1;
%         isi2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
%     end
%   
%     figure;hist(isi2,[1:states]);



%% seperate common spikes
for rr = 0
    B=BinningSpike(38,:);
    A=BinningSpike(6,:);
    delay = 10*rr;
    tolerence = 20;
    B1=zeros(1,length(B));
    B2=B;
    
    for i = 1:length(A)
        if A(i) ~= 0 
        B1(i-(delay+tolerence)/bin:i-(delay-tolerence)/bin)=B(i-(delay+tolerence)/bin:i-(delay-tolerence)/bin);
        B2(i-(delay+tolerence)/bin:i-(delay-tolerence)/bin)=0;
        else
        end
    end
    A1=zeros(1,length(A));
    A2=A;
    for i = 1:length(B)
        if B(i) ~= 0 
        A1(i-(delay+tolerence)/bin:i-(delay-tolerence)/bin)=A(i-(delay+tolerence)/bin:i-(delay-tolerence)/bin);
        A2(i-(delay+tolerence)/bin:i-(delay-tolerence)/bin)=0;
        else
        end
    end

%% Mutual Information
    for nn=1:5
        n=[A',B',B2',B1',A2'];
        Neurons = n(:,nn)';


%         n=roi(nn);Neurons = BinningSpike(n,:);

%         %% shuffle
%         r=randperm(length(Neurons));
%         for j=1:length(r)            
%             sNeurons(j)=Neurons(r(j));
%         end
%         Neurons=sNeurons;

    backward=ceil(1000/bin); forward=ceil(1000/bin);
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

%     informationp(temp) = ImExtrapolation_function(dat,8);
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
        c=corrcoef(x,y);
        corrf(temp)=c(2,1);
%       
%     informationf(temp) = ImExtrapolation_function(dat,8);
    end

    information=[informationp informationf]/BinningInterval;
    corr=[corrp corrf];
%     fr(z)=mean(Neurons)/BinningInterval; %firing rate(Hz)
    
    t=[-backward*bin:bin:forward*bin];  
    figure(14);hold on;plot(t,information,'LineWidth',2);%'color',cc(z,:),,'DisplayName',num2str(n)
    l={'N','P','P"','P1','N"'};legendInfo{nn} = [l{nn}];legend(legendInfo);
    set(gca,'FontSize',16)
    xlabel('\deltat (ms)');ylabel('MI (bits/s)');
    
        rrr(:,nn)=information';
    figure(15);hold on;plot(t,corr,'LineWidth',2);
    end
    
% pks(rr)=max(information);  

end
end
