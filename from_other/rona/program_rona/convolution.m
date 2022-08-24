% Convoulted real STA with original stimulaiton, by Rona, 2018
load('E:\google_rona\20170112\STA\STAcSTA_20170112_43.mat');window =0.8;window2 = 0;BinningInterval = 0.001;
% load('E:\google_rona\20170524\test4.mat');nSTA=STA';
% nSTA = -1.5*gampdf(0:1000,20,20)+gampdf(0:1000,12,12);
% nSTA = nSTA./max(abs(nSTA));
% figure; plot(nSTA)
% BinningInterval=10*10^-3;
SamplingRate=20000;
% 
test = repmat([nSTA],1,SamplingRate/1000);
r = reshape(test,[length(nSTA),SamplingRate/1000])';
NeuronsSTA=[];STA=[];
STA = r(:);
STA =STA/max(STA);
STA = STA-mean(STA);
NeuronsSTA = conv(fliplr(STA),isi);
NeuronsSTA = NeuronsSTA(SamplingRate/1000*window2/BinningInterval+1:end);
NeuronsSTA = NeuronsSTA/max(NeuronsSTA);
region = SamplingRate/1000*window/BinningInterval+500:length(isi)-SamplingRate/1000*window/BinningInterval-500;

isi2 = downsample(isi(region),0.01*SamplingRate);
NeuronsSTA2 =  downsample(NeuronsSTA(region),0.01*SamplingRate);
% isi2 = isi(region);
% NeuronsSTA2 = NeuronsSTA(region);

% NeuronsSTA = conv(inten,fliplr(STA));
% NeuronsSTA=NeuronsSTA(1:length(inten));
% NeuronsSTA=NeuronsSTA/max(NeuronsSTA);



                                                                                                                       
% for i=1:length(NeuronsSTA)
%     if NeuronsSTA(i)<mean(NeuronsSTA);
%         NeuronsSTA3(i)=1;
%     else
%         NeuronsSTA3(i)=0;
%     end
% end



% NeuronsSTA2=[];
% states=8;
% X=NeuronsSTA;
% nX = sort(X);
% abin = length(nX)/states;
% intervals = [nX(1:abin:end) inf]; %[-inf m-45 m-35 m-25 m-15 m-5 m+5 m+15 m+25 m+35 m+45 inf]; %[-inf min(X) adaptive max(X) inf];%[-inf min(isi):5:max(isi) inf];%
% temp=0;
% for jj = 1:BinningInterval*SamplingRate:length(X)
%     temp=temp+1;
% %         NeuronsSTA2(temp) = find(X(jj)<intervals,1)-1; % stimulus for every 50ms
%     NeuronsSTA2(temp) =  NeuronsSTA(jj);
% end
% figure;hist(NeuronsSTA2,states);

[t1,mi2]=MIfun(NeuronsSTA2,isi2);
% figure(111);hold on;plot(t1,mi1,'b',t1,mi2,'r');
figure(225);hold on;plot(t1,mi2)
