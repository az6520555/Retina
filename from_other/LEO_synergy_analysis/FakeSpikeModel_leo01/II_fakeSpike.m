clear all;
close all;
%% you can load different trajectory
load('merge_0224_HMM_UL_DR_G4.5_5min_Q100_6.5mW.mat')
%% Calculate correlation time. Need to install a toolbox. You can ignore this.
% acf = autocorr(bin_pos,100);
% corr_time = interp1(acf,1:length(acf),0.5,'linear')/60;
%% change dt here
dt = 0.5; %s
x = bin_pos(1:end);
v = finite_diff(x ,4);%v here is actually dx
r = x+v*dt*60; %%Generate spikes from poisson process
RespoSN = 4;
thesholds = [[std(r) 2*std(r) 3*std(r)]+mean(r) max(r)];
isir = [];
for jj=1:length(r)
    isir(jj) = find(r(jj)<=thesholds,1);
end
% 
% r = v;
% thesholds = [[std(r) 2*std(r) 3*std(r)]+mean(r) max(r)];
% for jj=1:length(r)
%     isir(jj) = isir(jj) + find(r(jj)<=thesholds,1);
% end
% r = x;
% thesholds = [[std(r) 2*std(r) 3*std(r)]+mean(r) max(r)];
% for jj=1:length(r)
%     isir(jj) = isir(jj) + find(r(jj)<=thesholds,1);
% end

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
legend('MI(x, �?��)','MI(v, �?��)', 'MI([x,v], �?��)', 'MI(x, �?��)+MI(v, �?��)');

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