% N P cell statistic
clear all
close all
path='G:\我的雲端硬碟\Retina exp\exp data\整理\NP_classification\OU\';
cd(path);
all_file = dir('*.mat');
P_cell=cell(1,length(all_file));
a=[1 2 3];b=[4 5 6];c=[7 8 9];d=[10 11];e=[12 13];f=[14 15];g=[16];
for z = 1:length(all_file)
    load([path,all_file(z).name])
    P_cell{z}=P_channel;
    N_cell{z}=N_channel;
end
for i=1:length(P_cell)
    P_num(i)=length(P_cell{i});
    N_num(i)=length(N_cell{i});
end
pa=sum(P_num(a));pb=sum(P_num(b));pc=sum(P_num(c));pd=sum(P_num(d));
pe=sum(P_num(e));pf=sum(P_num(f));pg=sum(P_num(g));
na=sum(N_num(a));nb=sum(N_num(b));nc=sum(N_num(c));nd=sum(N_num(d));
ne=sum(N_num(e));nf=sum(N_num(f));ng=sum(N_num(g));
ppp=[pa pb pc pd pe pf pg];
nnn=[na nb nc nd ne nf ng];
histogram('BinEdges',0:7,'BinCounts',ppp);hold on
histogram('BinEdges',0:7,'BinCounts',nnn)

P_tot=sum(P_num);
N_tot=sum(N_num);
tot_cell=P_tot+N_tot;

