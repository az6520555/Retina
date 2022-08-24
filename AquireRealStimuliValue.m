% transfer a_data to actual output voltage
clear
cd('E:\Chou\20190103\')
% load('20200116_cSTA_dt=50ms.mat')
filename='20180103_HMM_G=20.mcd';
file = strcat(filename);
AllDataInfo =datastrm(file) ;

MicrovoltsPerAD = getfield(AllDataInfo,'MicrovoltsPerAD2') %microvolt per digit
ZeroADValue = getfield(AllDataInfo,'ZeroADValue') % zero value in digit

sti=a_data(1,:); % original stimuli
sti=(sti-ZeroADValue).*MicrovoltsPerAD(2)*10^(-6);
