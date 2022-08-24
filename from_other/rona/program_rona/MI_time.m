% clearvars -except m sp1 sp2 isi2 statistime statispks STA w ro
clear all
close all
cd('\\192.168.0.100\Experiment\Retina\Rona\Exp\20180320\retina2\tHMM');

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
for z = [1:3]
%     clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA ro
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    bin=50;  BinningInterval = bin*10^-3; 
    backward=ceil(3000/bin); forward=ceil(3000/bin);
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
%     if length(TimeStamps)==1
%         TimeStamps(2)=TimeStamps(1)+290;
%     end
%     
   %% diode as isi2
%    load(['E:\google_rona\20171023\diode\diode_',filename]);
%      [~,locs_a2]=findpeaks(diff(diff(a2)),'MINPEAKHEIGHT',5*std(diff(diff(a2))));
%      TimeStamps_a2 = (locs_a2)/SamplingRate;
%      
% %     [b,a]=butter(2,60/20000,'low');
% %     a_data2=filter(b,a,callumin)';
%     a_data2=eyf;
%     isi = a_data2(TimeStamps_a2(1)*20000:TimeStamps_a2(end)*20000);
% %     figure(z);autocorr(isi,100000);
     %% Spike process
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end
    [n,out] = hist(TimeStamps,BinningTime);
    Stimuli = n;
    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;% figure;plot(BinningTime,sum(BinningSpike),BinningTime,10*Stimuli,'o')
%     figure;imagesc(BinningTime,[1:60],BinningSpike)
    
  
%% timewindow of stim %%%
    isi = TimeStamps(2:1:end) - TimeStamps(1:1:end-1); isi = isi*1000;  %for pseudo-period
    X = zeros(size(BinningTime));
    temp = 1;
    for ii = 1:length(X)
        if BinningTime(ii)<=TimeStamps(temp)
            X(ii) = isi(temp);%find(isi<TimeStamps(temp),1);
        else
            temp = temp + 1;
            if temp>length(isi); break; end
            X(ii) = isi(temp);
        end
    end

    m = mean(isi);  %% tt = abs(isi-m);  dev = max(tt(find(tt<300)));
    isi2=[];
    states=25;
    X(find(X==0))=mean(X);
    nX = sort(X);
    abin = length(nX)/states;
    intervals = [nX(1:abin:end) inf]; 
    for jj = 1:length(X)
        [a,b] = find(X(jj)<intervals,1); % stimulus for every 50ms
        isi2(jj) = b-1;
    end
  %% state of changing rate %%% 
%     isi3=[];isi4=[];
%     temp=0;
%     for j=1:length(isi)
%         temp=temp+1;
%         isi3(temp) = isi(j);
%     end
%     
%     for i=1:length(isi3)-1
%         isi4(i)=isi3(i+1)-isi3(i);
%     end
%     
%     X = zeros(size(BinningTime));
%     temp = 1;
%     for ii = 1:length(X)
%         if BinningTime(ii)<=TimeStamps(temp)
%             X(ii) = isi4(temp);%find(isi<TimeStamps(temp),1);
%         else
%             temp = temp + 1;
%             if temp>length(isi4); break; end
%             X(ii) = isi4(temp);
%         end
%     end
%     
%     isi2=[];
%     states=25;
%     nX = sort(X);
%     abin = length(nX)/states;
%     intervals = [nX(1:abin:end) inf]; 
%     temp=0;
%     for jj = 1:length(X)
%         temp=temp+1;
%         isi2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
%     end 
%       figure;hist(isi2,[1:states]);

%% Mutual Information
    for nn=1:length(roi)
        n=roi(nn);
        Neurons = BinningSpike(n,:); 
         %% shuffle
%         r=randperm(length(Neurons));
%         for j=1:length(r)            
%             sNeurons(j)=Neurons(r(j));
%         end
%         Neurons=sNeurons;

    dat=[];informationp=[];temp=backward+2;
    for i=1:backward+1 %past(t<0)

        if length(Neurons)==length(isi2)
            x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
        else
            x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1)-1)';
        end    
        y=isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        [N,C]=hist3(dat{i},[max(Neurons),25]); %20:dividing firing rate  6:# of stim
        px=sum(N,1)/sum(sum(N)); % x:stim
        py=sum(N,2)/sum(sum(N)); % y:word
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
              temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
            end
        end
        temp=temp-1;
        informationp(temp)=nansum(temp2(:));
        c=corrcoef(x,y);
        corrp(temp)=c(2,1);
    end  

    dat=[];informationf=[];temp=0;sdat=[];
    for i=1:forward
        if length(Neurons)==length(isi2)
            x =Neurons(forward+1-i:length(Neurons)-backward-i)';
        else
            x = Neurons(forward+1-i:length(Neurons)-backward-i-1)';
        end     
        y = isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        [N,C]=hist3(dat{i},[max(Neurons),25]); %20:dividing firing rate  6:# of stim
        px=sum(N,1)/sum(sum(N)); % x:stim
        py=sum(N,2)/sum(sum(N)); % y:word
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
                temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
            end
        end
        temp=temp+1;
        informationf(temp)=nansum(temp2(:)); 
        c=corrcoef(x,y);
        corrf(temp)=c(2,1);
    end

    information=[informationp informationf]/BinningInterval;
    corr=[corrp corrf];
%     fr(z)=mean(Neurons)/BinningInterval; %firing rate(Hz)
    [pks(z,nn),plocs(z,nn)]=max(information);
    t=[-backward*bin:bin:forward*bin];  
    figure(100);hold on;subplot(8,8,rr(n));
    xlim([-backward*bin,forward*bin]);
%     figure(n);hold on;
    plot(t,information,'LineWidth',2,'LineStyle','-');%,'color',cc(z,:)

%     set(gca,'XTick',[])
    %     l={'A','B','B"','B1'};legendInfo{nn} = [l{nn}];legend(legendInfo);
%     set(gca,'FontSize',16)
%     xlabel('\deltat (ms)');ylabel('MI (bits)');

    end
  end

pks(find(pks==0))=NaN;
% figure;errorbar(nanmean(pks'),nanstd(pks'),'--o');
% statistime = (-backward+plocs)*bin;
% figure;errorbar([50,300,600],fliplr(mean(statistime')),fliplr(std(statistime')),'--o');
% figure;errorbar(mean(w),std(w),'--o');
% 