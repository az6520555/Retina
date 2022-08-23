function [ey,a2,ss]=LED_pattern(stimu,file,pw)

rate=20000;
%% new ONOFF
if stimu == 'oo'
    repeat = 40; 
    on_t = 0.5; %s
    off_t = 0.5; %s
    mean_t = 1.5; %s
    period = (2*mean_t+on_t+off_t);
    mean_lumin=10;
    mean_i = mean_lumin;
    i_offset =0.4*mean_i;
    at = 60;%adaptation time
    ey = zeros; 
    a2 = -0.1*ones;
    ey(1:at*rate)=mean_i;
    for j = (0:repeat-1) % how many a2 in a trial
        ey( (at+period*j)*rate+1 : (at+period*j+on_t)*rate ) = mean_i+i_offset;
        ey( (at+period*j+on_t)*rate+1 : (at+period*j+mean_t+on_t)*rate ) = mean_i;
        ey( (at+period*j+mean_t+on_t)*rate+1 : (at+period*j+mean_t+on_t+off_t)*rate ) = mean_i-i_offset;
        ey( (at+period*j+mean_t+on_t+off_t)*rate+1 : (at+period*j+2*mean_t+on_t+off_t)*rate ) = mean_i;
        a2( (at+period*j)*rate+1 : (at+period*j+on_t)*rate ) = 1;
        a2( (at+period*j+on_t)*rate+1 : (at+period*j+mean_t+on_t)*rate ) = 0;
        a2( (at+period*j+mean_t+on_t)*rate+1 : (at+period*j+mean_t+on_t+off_t)*rate ) = 1;
        a2( (at+period*j+mean_t+on_t+off_t)*rate+1 : (at+period*j+2*mean_t+on_t+off_t)*rate ) = 0;
    end
    x1 = 1:length(ey);
    figure;plot(x1,ey);hold on; plot(x1,a2)
    
    ss = ['_Gonoff_'];
   

%% ConatrastSTA for normal distribution
elseif stimu == 'cn'
    T = 300; %nth in sec
    unit = 0.05;  %time step unit in ~ms
    steps = T/unit;
    mean_lumin=10;
    d=mean_lumin*(randn(1,steps));
    d=d-min(d);
    d=d+mean(d)/5;
    d=d*(mean_lumin/mean(d));
    at=60;
    ey1=[];ey=[];ey0=[];
    for i=1:steps
        ey1(rate*unit*(i-1)+1:rate*unit*i)=d(i);
    end
    ey0=mean_lumin*ones(1,at*rate);%REST
    ey=[ey0 ey1]; 
    a2=[zeros(1,at*rate) ones(1,1*rate) zeros(1,length(ey)-(at+1)*rate)];
    a2((T+at-1)*rate:(T+at)*rate)=1;
    figure;plot(ey);
    ss = ['_cSTA_'];
    
%% OSR %%%%%%%
elseif stimu == 'os'
    rep = 20; %20 trials
    rest = 5; %5s inter-trial
    trial = 30; %20s trial
    ad=30;
    w = 0.05; %pulse width =50ms
    p = 0.20; %period 
    allX(1:2:2*6000)= w;
    allX(2:2:2*6000)= p-w;%X; for OSR!
    ey1=[];
    for i=1:rep*2 %length(allX) for OSR!
        if mod(i,2)==0
           ey1=[ey1 zeros(1,round(allX(i)*rate))];
        else
           ey1=[ey1 ones(1,round(allX(i)*rate))];
        end
    end
    ey = [];
    for i = 1:trial
        ey = [ey ey1 zeros(1,round(rest*rate))];
    end
   
    a2 = ey;
    ey = 10-10*ey;
    ey = [10*ones(1,ad*rate) ey];
    a2 = [10*zeros(1,ad*rate) a2];
    figure; plot(ey);
    length(ey)/20000
    ss = ['_OSR_',num2str(p*1000),'_'];
%% OSR sudden change frequency recording different changing frequencies
elseif stimu == 'oc'
    rep = 20; %20 trials
    rest = 5; %5s inter-trial
    trial = 20; %20s trial
    pw = 0.05; %pulse width =50ms
    p0 = 0.20; %period 
    dp = 0.01; %value of pulse frequency changed every time
    A = [ones(1,round(pw*rate)) zeros(1,round((p0-pw)*rate))]; %one pulse for original frequency
    for k=1:trial
        B{k} = [ones(1,round(pw*rate)) zeros(1,round(((p0-pw)+dp*(k-1))*rate))]; %one pulse for increased frequency
    end
    
    ey1=[];  
  
    for k=1:trial
        ey1=[];
        for j=1:5
            for i=1:rep
                ey1=[ey1 A];
            end
            for i=1:rep
                ey1=[ey1 B{k}];
            end
            ey1=[ey1 ones(1,rest*rate)];
        end
        singleperiod{k}=ey1 ;
    end
    ey=[];
    for i=1:trial
    ey=[ey singleperiod{i}];
    end
   
    a2 = ey;
    ey = 0.15*ey;
    figure; plot(ey);
    length(ey)/20000
    ss = ['_OSR_SuddenChangeFrequency_',[num2str(p0*1000),'to',num2str((p0+dp*(trial-1))*1000)],'_'];
%% OSR sudden change frequency for only one second frequency
elseif stimu == 'ol'
    rep = 20; %20 pulses
    rest = 5; %5s inter-trial
    trial = 15; %20s trial
    pw = 0.05; %pulse width =50ms
    p0 = 0.20; %period 
    dp = 0.10; %value of pulse frequency changed every time
    tic
    A = [ones(1,round(pw*rate)) zeros(1,round((p0-pw)*rate))]; %one pulse for original frequency
    B = [ones(1,round(pw*rate)) zeros(1,round(((p0-pw)+dp)*rate))]; %one pulse for increased frequency
    ey1=[];
    
    for i=1:rep
        ey1=[ey1 A];
    end
    for i=1:rep
        ey1=[ey1 B];
    end
    ey1=[ey1 zeros(1,rest*rate)];
    ey=[];
    for i=1:trial
        ey=[ey ey1];
    end
    toc
    
    a2 = ey;
    ey = 0.15*ey;
    figure; plot(ey);
    length(ey)/20000
    ss = ['_OSR_SuddenChangeFrequency_',[num2str(p0*1000),'to',num2str((p0+dp)*1000)],'_'];
%% Continuous increasing frequency
elseif stimu == 'jk'
    fl=5;
    fh=30;
    duty=1/2;  % might be changeable
    num_pulse=5000;
    df=(fh-fl)/num_pulse;
    I=10;
    beforetime=60;
    restperiod=120;
    repeat=2;
    f=fl;
    t=[];
    beforetest=I*ones(1,round(beforetime*rate));
    rest=zeros(1,round(restperiod*rate));
    a2before=zeros(1,length(beforetest));
    a2rest=zeros(1,length(rest));
    sti=[];
    for i=1:num_pulse
        f=f+df;
        p=1/f;
    	sti=[sti I*1.5*ones(1,round(p*duty*rate)) I*0.5*ones(1,round(p*(1-duty)*rate))];
    end
    aa2=zeros(1,length(sti));
    aa2(1:rate)=1;
    aa2(end-rate:end)=1;
    
    ey=[];
    a2=[];
    for j=1:repeat-1
        ey=[ey beforetest sti rest];
        a2=[a2 a2before aa2 a2rest];
    end
    ey=[ey beforetest sti];
    a2=[a2 a2before aa2];
    
    t=1/rate:1/rate:length(ey)/rate;
    plot(t,ey)

    ss=['_DecPeriodPulse_',num2str(fl),'hz_to',num2str(fh),'hz_duty',num2str(duty)];
    

%% constant period pulse (onoff stimulus)
elseif stimu == 'cp'
    f=8;
    p=1/f;
    pw=1/2*p;
    n=500;
    C=1;
    repeat=3;
    adaptation=60;
    rest=60;
    I_high=0.2;
    I_low=(1-C)/(1+C)*I_high;
    A=ones(1,round(pw*rate));
    B=ones(1,round((p-pw)*rate));
    ey1=[];
    for i=1:n
        ey1=[ey1 I_high*A I_low*B];
    end
    ey=[];
    ey=[ey 0.5*I_high*ones(1,round(adaptation*rate))];
    for i=1:repeat-1
        ey=[ey ey1 I_low*ones(1,rest*rate) 0.5*I_high*ones(1,adaptation*rate)];
    end
    ey=[ey ey1];
    
    a2=ey;
    t=1/rate:1/rate:length(ey)/rate;
    plot(t,ey)
    ss=['ConstantPulseOn_f',num2str(f),'_C',num2str(C),'_']
     
%% constant period pulse (constant mean light intensity)
elseif stimu == 'lp'    
    sti_t=30;
    f=9;
    p=1/f;
    pw=0.05;
    n=fix(sti_t/p);
    repeat=5;
    rest=20;
    C=0.3;
    I_mean=0.1;
    I_high=(1+C)*I_mean;
    I_low=(1-C)*I_mean;
    A=ones(1,round(pw*rate));
    B=ones(1,round((p-pw)*rate));
    ey1=[];
    for i=1:n
        ey1=[ey1 I_high*A I_low*B];
    end
    ey=[];
    for i=1:repeat-1
        ey=[ey ey1 I_low*ones(1,rest*rate)];
    end
    ey=[ey ey1];
    
    a2=ey;
    t=1/rate:1/rate:length(ey)/rate;
    plot(t,ey)
    ss=['ConstantPulse_f',num2str(f),'_C',num2str(C),'_']
    
%% constant period pulse (off stimulus)
elseif stimu == 'cf' 
    f=10;
    p=1/f;
    pw=p/2;
    n=32;
    C=1;
    adapt_time=30;
    repeat=50;
    rest=2;
    I_high=0.2;
    I_low=(1-C)/(1+C)*I_high;
    A=ones(1,round(pw*rate));
    B=ones(1,round((p-pw)*rate));
    ey1=[];
    for i=1:n
        ey1=[ey1 I_low*A I_high*B];
    end
    ey=[I_high*ones(1,adapt_time*rate)];
    for i=1:repeat-1
        ey=[ey ey1 I_high*ones(1,rest*rate)];
    end
    ey=[ey ey1];
    
    a2=-ey+I_high;
    t=1/rate:1/rate:length(ey)/rate;
    plot(t,ey)
    ss=['ConstantPulseOn_f',num2str(f),'_C',num2str(C),'_']
    
%% Load mat file
elseif stimu == 'ld'
%     path='\\192.168.1.100\Experiment\Retina\Chou\stimulation code\20200318_OU_filtering\';
%     filename='OU_original_tau=0p5.mat'
    load(file);
    
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    try 
        ey=inten;
    catch
    end

%     figure(3);plot(t,ey);hold on 
    ss=[name,'_'];
%% Gaussian pulse (ON)
elseif stimu == 'gn'
    repeat=41;
    rest=2;
    at=60;
    a2du=0.5;
    minI=10; % mW/m^2
    amp=minI*0.4;
    xp=-pw*5*rate:pw*5*rate;
    ey1=[];
    ey1=amp*exp(-(xp/(rate*sqrt(2)*(pw/2.355))).^2);
    ey1=ey1+minI;
    ey=[];
    a2=[];
    ey=minI*ones(1,rate*at);
    a2=zeros(1,length(ey));
    for i=1:repeat-1
        ey=[ey ey1 minI*ones(1,rest*rate)];
        a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)+(rest-a2du)*rate)];
    end
    ey=[ey ey1];
    a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)-a2du*rate)];
    
    ss=['_GaussianPulseOn_pulsewidth=',num2str(pw),'_'];
    
%% Gaussian pulse
elseif stimu == 'gf'
%     pw=1;%
    repeat=41;
    rest=2;
    at=60;
    a2du=0.5;
    maxI=10;
    amp=maxI*0.4;
    xp=-pw*5*rate:pw*5*rate;
    ey1=[];
    ey1=maxI-amp*exp(-(xp/(rate*sqrt(2)*(pw/2.355))).^2);
    ey1=ey1;
    ey=[];
    a2=[];
    ey=maxI*ones(1,rate*at);
    a2=zeros(1,length(ey));
    for i=1:repeat-1
        ey=[ey ey1 maxI*ones(1,rest*rate)];
        a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)+(rest-a2du)*rate)];
    end
    ey=[ey ey1];
    a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)-a2du*rate)];
    
    ss=['_GaussianPulseOff_pulsewidth=',num2str(pw),'_'];

    %% Gaussian pulse ON sudden stop
elseif stimu=='nt'
    repeat=41;
    rest=5;
    at=60;
    a2du=0.5;
    maxI=10;
    amp=maxI*0.4;
    xp=-pw*5*rate:0;
    ey1=[];
    ey1=maxI+amp*exp(-(xp/(rate*sqrt(2)*(pw/2.355))).^2);
    ey=[];
    a2=[];
    ey=maxI*ones(1,rate*at);
    a2=zeros(1,length(ey));
    for i=1:repeat-1
        ey=[ey ey1 maxI*ones(1,rest*rate)];
        a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)+(rest-a2du)*rate)];
    end
    ey=[ey ey1];
    a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)-a2du*rate)];
    ss=['_GaussianON_SuddenStop_pulsewidth=',num2str(pw),'_'];
    
    %% Gaussian pulse ON sustain
elseif stimu=='ns'
    repeat=41;
    rest=5;
    sustain=2;
    at=60;
    a2du=0.5;
    maxI=10;
    amp=maxI*0.4;
    xp=-pw*5*rate:0;
    ey1=[];
    ey1=maxI+amp*exp(-(xp/(rate*sqrt(2)*(pw/2.355))).^2);
    ey1=[ey1 (maxI+amp)*ones(1,sustain*rate)];
    ey=[];
    a2=[];
    ey=maxI*ones(1,rate*at);
    a2=zeros(1,length(ey));
    for i=1:repeat-1
        ey=[ey ey1 maxI*ones(1,rest*rate)];
        a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)+(rest-a2du)*rate)];
    end
    ey=[ey ey1];
    a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)-a2du*rate)];
    
    ss=['_GaussianON_Sustain_pulsewidth=',num2str(pw),'_'];
    
    %% Gaussian pulse OFF sustain
elseif stimu=='fs'
    repeat=41;
    rest=5;
    sustain=2;
    at=60;
    a2du=0.5;
    maxI=10;
    amp=maxI*0.4;
    xp=-pw*5*rate:0;
    ey1=[];
    ey1=maxI-amp*exp(-(xp/(rate*sqrt(2)*(pw/2.355))).^2);
    ey1=[ey1 (maxI-amp)*ones(1,sustain*rate)];
    ey=[];
    a2=[];
    ey=maxI*ones(1,rate*at);
    a2=zeros(1,length(ey));
    for i=1:repeat-1
        ey=[ey ey1 maxI*ones(1,rest*rate)];
        a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)+(rest-a2du)*rate)];
    end
    ey=[ey ey1];
    a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)-a2du*rate)];
    
    ss=['_GaussianOFF_Sustain_pulsewidth=',num2str(pw),'_'];
    
    %% Gaussian pulse ON sudden stop
elseif stimu=='ft'
    repeat=41;
    rest=5;
    at=60;
    a2du=0.5;
    maxI=10;
    amp=maxI*0.4;
    xp=-pw*5*rate:0;
    ey1=[];
    ey1=maxI-amp*exp(-(xp/(rate*sqrt(2)*(pw/2.355))).^2);
    ey1=ey1;
    ey=[];
    a2=[];
    ey=maxI*ones(1,rate*at);
    a2=zeros(1,length(ey));
    for i=1:repeat-1
        ey=[ey ey1 maxI*ones(1,rest*rate)];
        a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)+(rest-a2du)*rate)];
    end
    ey=[ey ey1];
    a2=[a2 ones(1,a2du*rate) zeros(1,length(ey1)-a2du*rate)];
    
    ss=['_GaussianOFF_SuddenStop_pulsewidth=',num2str(pw),'_'];
    
    
end