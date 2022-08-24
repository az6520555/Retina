%% stimulus with different intensity pulse
clear all
close all
cd('E:\rona\20170330');
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
SamplingRate=20000;
cc=hsv(n_file);

for z =11:13
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA ro
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';

    bin=10;BinningInterval=bin*10^-3; 
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    backward=ceil(3000/bin); forward=ceil(3000/bin);
    roi=[35];
    
    a_data2 = a_data(3,:);   
    a_data2 = a_data2 - a_data2(1);
    SamplingRate=20000;
    [~,locs]=findpeaks(diff(a_data2),'MINPEAKHEIGHT',5*std(diff(a_data2)));
    analog_loc = (locs)/SamplingRate;
    TimeStamps = analog_loc;

    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end
    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;
    figure;imagesc(BinningTime,[1:60],BinningSpike)
    
    for nn=1:length(roi)
        Neurons = BinningSpike(roi(nn),:);
        a_data2=a_data(3,:);
        % isi = a_data2(TimeStamps(1)*20000+1000:BinningInterval*20000*6:TimeStamps(end)*20000+1000);% figure;plot(isi);
        for i=1:length(TimeStamps)
            isi(i) = a_data2(ceil(TimeStamps(i)*20000+200));
        end
        m = mean(isi); % mV
        DataTime = TimeStamps(end)-TimeStamps(1); %s
        % %     figure;autocorr(isi,20000);
        for i=1:length(TimeStamps)
              t(i)=TimeStamps(i)*20000;
        end
        figure;plot(a_data2); 
        hold on; plot(t,isi,'*');
     %% state of light intensity %%%
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
        isi2=[];
        states=8;
        %     X=isi;
        nX = sort(X);
        abin = length(nX)/states;
        intervals = [nX(1:abin:end) inf]; 
        temp=0;
        for jj = 1:length(X)
            temp=temp+1;
            isi2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
        end
        figure;hist(isi2);

    %% shuffle
    % r=randperm(length(Neurons));
    % for j=1:length(r)            
    %     sNeurons(j)=Neurons(r(j));
    % end
    % Neurons=sNeurons; 
    %%
    temp=0;
    dat=[];informationp=[];temp=backward+2;
    for i=1:backward+1 %past(t<0)
        x=Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
        y=isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
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
        x = Neurons(forward+1-i:length(Neurons)-backward-i)';
        y = isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
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
    information=[informationp informationf];
    corr=[corrp corrf];
    t=[-backward*bin:bin:forward*bin];  
    figure(roi(nn));hold on;plot(t,information,'LineWidth',2);%'color',cc(z,:),
    set(gca,'FontSize',16)
    xlabel('t-t (ms)');ylabel('MI (bits)');
    figure(10);hold on;plot(t,corr);
    end
end
    


    