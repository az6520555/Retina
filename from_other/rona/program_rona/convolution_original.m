load('E:\google_rona\20170112\STA\STAcSTA_20170112_43.mat')

% nSTA = -1.5*gampdf(0:1000,20,20)+gampdf(0:1000,12,12);
% nSTA = nSTA./max(abs(nSTA));
% figure; plot(nSTA)

BinningInterval=10*10^-3;
SamplingRate=20000;
test = repmat([nSTA],1,SamplingRate/1000);
r = reshape(test,[length(nSTA),SamplingRate/1000])';
NeuronsSTA=[];STA=[];
STA = r(:);
% STA =STA/max(STA);
NeuronsSTA = conv(fliplr(STA),isi);
% NeuronsSTA = conv(STA,isi);
NeuronsSTA=NeuronsSTA(1:length(isi));
NeuronsSTA=NeuronsSTA/max(NeuronsSTA);


% for i=1:length(NeuronsSTA)
%     if NeuronsSTA(i)<mean(NeuronsSTA);
%         NeuronsSTA3(i)=1;
%     else
%         NeuronsSTA3(i)=0;
%     end
% end


% NeuronsSTA=NeuronsSTA(1:100:end-19960);
% NeuronsSTA=NeuronsSTA/max(NeuronsSTA);
% NeuronsSTA=NeuronsSTA;

NeuronsSTA2=[];
states=8;
X=NeuronsSTA;
nX = sort(X);
abin = length(nX)/states;
intervals = [nX(1:abin:end) inf]; %[-inf m-45 m-35 m-25 m-15 m-5 m+5 m+15 m+25 m+35 m+45 inf]; %[-inf min(X) adaptive max(X) inf];%[-inf min(isi):5:max(isi) inf];%
temp=0;
for jj = 1:BinningInterval*SamplingRate:length(X)
    temp=temp+1;
%         NeuronsSTA2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
    NeuronsSTA2(temp) =  NeuronsSTA(jj);
end
figure;hist(NeuronsSTA2,states);

%% shuffle
%         r=randperm(length(Neurons));
%         for j=1:length(r)            
%             sNeurons(j)=Neurons(r(j));
%         end
%         NeuronsSTA=sNeurons;
bin = 10;  BinningInterval = bin*10^-3;  %ms
backward=ceil(1000/bin); forward=ceil(1000/bin);
dat=[];informationp=[];temp=backward+2;corrp=[];
    for i=1:backward+1 %past(t<0)
        x=NeuronsSTA2((i-1)+forward+1:length(NeuronsSTA2)-backward+(i-1))';
        y=isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        norm=1;

        [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
        px=sum(N,1)/sum(sum(N)); % x:stim
        py=sum(N,2)/sum(sum(N)); % y:word
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
              temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2)/norm;
            end
        end
        temp=temp-1;
        informationp(temp)=nansum(temp2(:));
        c=corrcoef(x,y);
        corrp(temp)=c(2,1);
    end  

    dat=[];informationf=[];temp=0;sdat=[];corrf=[];
    for i=1:forward
        x = NeuronsSTA2(forward+1-i:length(NeuronsSTA2)-backward-i)';
        y = isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        norm=1;

        [N,C]=hist3(dat{i}); %20:dividing firing rate  6:# of stim
        px=sum(N,1)/sum(sum(N)); % x:stim
        py=sum(N,2)/sum(sum(N)); % y:word
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
                temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2)/norm;
            end
        end
        temp=temp+1;
        informationf(temp)=nansum(temp2(:)); 
        c=corrcoef(x,y);
        corrf(temp)=c(2,1);
    end
    
    information=[informationp informationf]/BinningInterval;
    corr=[];
    corr=[corrp corrf];

    t=[-backward*bin:bin:forward*bin];
    name(name=='_')=' ';
    figure(n);hold on;plot(t,information,'color',cc(z,:),'LineWidth',2);hold on
    xlabel('time (ms)');ylabel('MI (bits/s)');
    title([name,' n=',num2str(n)]);
    
%     figure(8);plot(t,corr,'color',cc(z,:),'LineWidth',2);hold on
%     xlabel('time (ms)');ylabel('cross corr');
%     title([name,' n=',num2str(n)]);