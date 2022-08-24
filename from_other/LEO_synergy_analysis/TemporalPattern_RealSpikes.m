clear all;
close all;
%% you can load different trajectory 
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\analysis code\LEO_synergy_analysis')
load('20200408_OU_cutoff=4_sort_unit2.mat')
%% Calculate correlation time. Need to install a toolbox. You can ignore this.
% acf = autocorr(bin_pos,100);
% corr_time = interp1(acf,1:length(acf),0.5,'linear')/60;
%% change dt here
BinningInterval=0.01;

[b,a] = butter(2,50/20000,'low'); % set butter filter
a_data2 = filter(b,a,a_data(1,:));
x = a_data2(TimeStamps(1)*20000:TimeStamps(length(TimeStamps))*20000);% figure;plot(isi);
x = downsample(x,200);
v = finite_diff(x ,4);%v here is actually dx

BinningTime = [TimeStamps(1) : BinningInterval : TimeStamps(end)];
[n,xout] = hist(Spikes{36},BinningTime);
isir=n+1;isir(1)=0;isir(end)=0;

StimuSN = 6;
sqrtStimuSN = 6;

nX=sort(x);
abin=length(nX)/sqrtStimuSN;
intervals=[nX(abin:abin:end) inf]; % inf: the last term: for all rested values
temp=0; isix=[];
for jj=1:length(x)
    isix(jj) = find(x(jj)<=intervals,1);
end
nX=sort(v);
abin=length(nX)/sqrtStimuSN;
intervals=[nX(abin:abin:end) inf]; % inf: the last term: for all rested values
isiv=[];
for jj=1:length(v)
    isiv(jj) = find(v(jj)<=intervals,1);
end
isii = sqrtStimuSN*(isiv-1) + isix;


nX=sort(x);
abin=length(nX)/StimuSN;
intervals=[nX(abin:abin:end) inf]; % inf: the last term: for all rested values
temp=0; isix=[];
for jj=1:length(x)
    isix(jj) = find(x(jj)<=intervals,1);
end
nX=sort(v);
abin=length(nX)/StimuSN;
intervals=[nX(abin:abin:end) inf]; % inf: the last term: for all rested values
isiv=[];
for jj=1:length(v)
    isiv(jj) = find(v(jj)<=intervals,1);
end


bin = BinningInterval*1000;
backward=ceil(5000/bin);
forward=ceil(5000/bin);
time=[-backward*bin:bin:forward*bin];
[information_x_r information_v_r redundant_I] = PIfunc(isir,isix, isiv,BinningInterval,backward,forward);
information_i_r = MIfunc(isir,isii,BinningInterval,backward,forward);
figure(4);
plot(time,information_x_r, 'r');hold on;
plot(time,information_v_r, 'b')
plot(time,information_i_r, 'k')
plot(time,information_x_r+ information_v_r, 'm')
legend('MI(x, ??›¾)','MI(v, ??›¾)', 'MI([x,v], ??›¾)', 'MI(x, ??›¾)+MI(v, ??›¾)');

%
% computation of synergy here
%
z = information_x_r + information_v_r -information_i_r;
synergy_I = redundant_I-z;
PI_x = information_x_r - redundant_I;
PI_v = information_v_r - redundant_I;
figure(6);
plot(time,PI_x, 'r');hold on;
plot(time,PI_v, 'b')
plot(time,synergy_I, 'k')
plot(time,redundant_I, 'g')
legend('PI_x','PI_v', 'synergy I', 'redundant I');