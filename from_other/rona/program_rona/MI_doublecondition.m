% Double condition mutual information, eg. MI(r;(x,v))
clear all
% close all

cd('E:\google_rona\20170524\HMM') ; % the folder of the files
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
SamplingRate=20000;
cc=hsv(n_file);

for z = 3
    clearvars -except all_file n_file z SamplingRate cc ey isi2
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    % user's parameter
    bin =10;  BinningInterval = bin*10^-3;  %ms
    backward=ceil(500/bin); forward=ceil(500/bin);
    roi=5;

    a_data2=a_data(3,:);
    
    isi = a_data2(TimeStamps(1)*20000:TimeStamps(length(TimeStamps))*20000);% figure;plot(isi);
    m = mean(isi); % mV
    DataTime =TimeStamps(end)-TimeStamps(1); %s

%% Spike process
    BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
    BinningSpike = zeros(60,length(BinningTime));
    for i = 1:60
        [n,xout] = hist(Spikes{i},BinningTime) ;
        BinningSpike(i,:) = n ;
    end

    BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;% figure;plot(BinningTime,sum(BinningSpike),BinningTime,10*Stimuli,'o')
%     figure;imagesc(BinningTime,[1:60],BinningSpike)
   
%     Neurons = BinningSpike(roi,:);
    N = BinningSpike(3,:);
    P = BinningSpike(16,:);
    
    %% state of light intensity %%%
    I=[];
    states=8;
    X = isi;
    nX = sort(X);
    abin = length(nX)/states;
    intervals = [nX(1:abin:end) inf]; 
    temp=0;
    for j = 1:BinningInterval*SamplingRate:length(X)
        temp=temp+1;
        I(temp) = isi(j);
        sI(temp) = find(X(j)<intervals,1)-1; % stimulus for every 50ms
    end
%     figure;hist(sI,[1:states]);

 %% state of changing rate %%% 
%     dI=[];
%     temp=0;    
%     for j=1:length(I)-1
%         dI(j)=I(j+1)-I(j);
%     end
%         
%     sdI=[];
%     states=8;
%     X=dI;
%     nX = sort(X);
%     abin = length(nX)/states;
%     intervals = [nX(1:abin:end) inf]; 
%     temp=0;
%     for jj = 1:length(X)
%         temp=temp+1;
%         sdI(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
%     end
%     figure;hist(sdI,[1:states]);

     %% shuffle
%         r=randperm(length(Neurons));
%         for j=1:length(r)            
%             sNeurons(j)=Neurons(r(j));
%         end
%         Neurons=sNeurons;
%% double conditioned Mutual Information MI(resp;(cond1,cond2))
    informationp=[];temp3=backward+2;
    cond1 = N+1;
    cond2 = N+1;
    resp = sI;
    for shiftp=1:backward+1; %shiftp=1,alligned
        x = cond1((shiftp-1)+forward+1:length(sI)-backward+(shiftp-1))';
        y = cond2((shiftp-1)+forward+1:length(sI)-backward+(shiftp-1))';
        z = resp(forward+1:length(N)-backward)';
        test=[x,y,z];

        dat={};
        for i=1:max(x)
            for j=1:max(y)
                px=[];py=[];
                px=find(test(:,1)==i);
                temp=0;
                for k=1:length(px)
                    if test(px(k),2)==j;
                        temp=temp+1;
                        py(temp)=px(k);
                    end
                    dat{i,j}=test(py,3);
                end
            end
        end

        temp1=0;label=0;dat2=[];
        for i=1:size(dat,1)
            for j=1:size(dat,2)
                label=label+1;
                for k=1:length(dat{i,j})
                    temp1=temp1+1;
                    dat2(temp1,1)=label;
                    dat2(temp1,2)=dat{i,j}(k);
                end
            end
        end

        [num,pos]=hist3(dat2,[64,8]); %20:dividing firing rate  6:# of stim
        px=sum(num,1)/sum(sum(num)); % x:stim
        py=sum(num,2)/sum(sum(num)); % y:word
        pxy=num/sum(sum(num));
        
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
              temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
            end
        end
        temp3=temp3-1;
        informationp(temp3)=nansum(temp2(:));       
    end
    
    informationf=[];temp3=0;test=[];
    for shiftf=1:forward;
        x=cond1(forward+1-shiftf:length(sI)-backward-shiftf)';
        y=cond2(forward+1-shiftf:length(sI)-backward-shiftf)';
        test=[x,y,z];

        fr=[];isi2=[];temp=0;dat={};
        for i=1:max(x)
            for j=1:max(y)
                px=[];py=[];
                px=find(test(:,1)==i);
                temp=0;
                for k=1:length(px)
                    if test(px(k),2)==j;
                        temp=temp+1;
                        py(temp)=px(k);
                    end
                    dat{i,j}=test(py,3);
                    fr(i,j)=sum(test(py,3));
                end
            end
        end
%         figure;imagesc(fr);

        temp1=0;label=0;dat2=[];
        for i=1:size(dat,1)
            for j=1:size(dat,2)
                label=label+1;
                for k=1:length(dat{i,j})
                    temp1=temp1+1;
                    dat2(temp1,1)=label;
                    dat2(temp1,2)=dat{i,j}(k);
                end
            end
        end

        [num,pos]=hist3(dat2,[64,8]); %20:dividing firing rate  6:# of stim
        px=sum(num,1)/sum(sum(num)); % x:stim
        py=sum(num,2)/sum(sum(num)); % y:word
        pxy=num/sum(sum(num));
        
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
              temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);
            end
        end
        temp3=temp3+1;
        informationf(temp3)=nansum(temp2(:));
        
    end
        information=[informationp informationf]/BinningInterval;
        t=[-backward*bin:bin:forward*bin];
        figure(roi);hold on;plot(t,information);

end