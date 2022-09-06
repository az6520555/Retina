clear all;

ddd=date;
% diode_path_cal=['\\192.168.1.100\Experiment\Retina\Chou\stimulus saving\',ddd,'\calibration'];
% diode_path_cal=['\\192.168.1.100\Public\chou\stimulus_saving\',ddd,'\calibration'];
diode_path_cal=['\\192.168.1.102\Public\Retina\Chou\stimulus_saving',ddd,'\calibration'];

rate=10000;  
du=10;
delay=0.02;
volt = linspace(0,12,du*rate);
% ======= Important: DAQ read in data delays the output for 0.02s =======
volt = [zeros(1,delay*rate) volt];

daq_out = daqmx_Task('chan','Dev1/ao1' ,'rate',rate, 'Mode', 'f');
daq_out.write(0);
%%%%%%%%%%%% step1:background lumin %%%%%%%%%%%
% daq_out = daqmx_Task('chan','Dev1/ao1' ,'rate',rate, 'Mode', 'f');
% daq_out.write(0);
% daq = daqmx_Task('Dev1/ai0');
% bk = daq.read;


%%%%%%%%%%%% step2:calibration line %%%%%%%%%%%
% daq_out = daqmx_Task('chan','Dev1/ao1' ,'rate',rate, 'Mode', 'f');
daq_in = daqmx_Task('chan','Dev1/ai0' ,'rate',rate, 'SampleNum', du*rate);
% daq_in = daqmx_Task('chan','Dev1/ai0' ,'rate',rate, 'SampleNum', length(ey));

daq_out.write(volt); 
daq_in.read;    

daq_out = daqmx_Task('chan','Dev1/ao1' ,'rate',rate, 'Mode', 'f');
daq_out.write(0);
% daq_in.wait;
lumin=daq_in.data;
% [b,a]=butter(2,20/rate,'low');
% lumin_filter=filter(b,a,lumin);
span=500;
lumin_filter=smooth(lumin,span); % smooth looks better
lumin_filter(1:span)=mean(lumin_filter(span+1:2*span));

% ================ transfer voltage signal to light intensity ===========
offset=min(lumin_filter)
lumin_filter=lumin_filter-offset;
Ip=lumin_filter/10.421/10^6;
r=0.37;
P=Ip/r;
A=13*10^-6;
inten=P/A*1000; % unit: mW/m^2
lumin_filter=inten;
% =======================================================================
%  plot(volt);title(['voltage v.s. time']);
figure(1);plot(volt(delay*rate+1:end),lumin);title(['voltage  v.s. time']);
figure(2);hold on; plot(volt(delay*rate+1:end),lumin_filter);title(['filtered lumin (mW/m^2) v.s. time']);
% figure(4) ; plot(volt,lumin_filter);title(['voltage v.s. lumin']);

mkdir(diode_path_cal)
save([diode_path_cal,'\calibration_PAC_',ddd],'lumin','lumin_filter','rate','volt','offset')
