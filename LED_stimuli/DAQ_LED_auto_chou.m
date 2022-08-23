clearvars -except rate volt lumin_filter offset stimu pw path file path parameters stimu_set...
     file_set pw_set
daq_out = daqmx_Task('chan','Dev1/ao1' ,'rate',rate, 'Mode', 'f');
daq_out.write(0);
rate=20000;
% load('E:\rona\20160318\map');
%%%%%%%%%%%%% step3:calibration %%%%%%%%%%%%%%%%%%%
x=volt;
y=(lumin_filter)';

[ey,a2,ss]=LED_pattern_chou(stimu,file,pw);

eyf=ey;

[y, index] = unique(y);
ex = interp1(y,x(index),eyf,'linear');% ex=calibrate voltage
daq = daqmx_Task('chan','Dev1/ao1' ,'rate',rate, 'SampleNum', length(eyf));
% daq.stop
daq_out = daqmx_Task('chan','Dev1/ao0:1' ,'rate',rate, 'Mode', 'f');
daq_in = daqmx_Task('chan','Dev1/ai0' ,'rate',rate, 'SampleNum', length(eyf));

A=[];
A(:,1) = a2(1:length(eyf));
A(:,2) = ex(1:length(eyf));

%%%% output %%%%%%%%%%
daq_out.write(A); 
daq_in.read;
callumin=daq_in.data;
[b,a]=butter(2,10/rate,'low');
callumin_filter=filter(b,a,daq_in.data);
daq_in.wait;


close all;
daq_out = daqmx_Task('chan','Dev1/ao1' ,'rate',rate, 'Mode', 'f');
daq_out.write(0);

% ================ transfer voltage signal to light intensity ===========
callumin=callumin-offset;
Ip=callumin/10.421/10^6;
r=0.37;
P=Ip/r;
A=13*10^-6;
inten=P/A*1000; % unit: mW/m^2  
callumin=inten;
[b,a]=butter(2,2000/rate,'low')
callumin=filter(b,a,callumin);
% =======================================================================

t=1/rate:1/rate:length(ex)/rate;
figure(5) ; plot(t,ex);title(['calibrated voltage v.s. time']);
figure(6) ; plot(t,callumin);title(['lumin v.s. time']);hold on;plot(t-0.08,ey,'r');plot(t-0.08,eyf,'g');
% figure(7) ; plot(callumin_filter);title(['lumin v.s. time']);hold on;plot(ey,'r');plot(eyf,'g');

ddd=date;
diode_path=['\\192.168.1.100\Public\chou\stimulus_saving\',ddd,'\']
mkdir(diode_path)
cd(diode_path)
filesaving=['diode_',ss,ddd];
i=0;
n=87;
while n~=0
    n=exist([filesaving,'_',num2str(i),'.mat'],'file')
    if n==0
        save([filesaving,'_',num2str(i),'.mat'],'callumin','callumin_filter','a2','ex','ey','eyf');
    end
    i=i+1;
end