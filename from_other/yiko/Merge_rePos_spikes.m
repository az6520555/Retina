%load mat file first
clear all
close all
%load ideal stimulation array
% cd('D:\Yiko\Files for Thesis\04122018') ; % for saving
G=3;
Dirti='RL';

%% load bar position worksapce
% load('D:\Yiko\HMM video frame\workspace_HMM RL G5 5min Q85')
% load('\\192.168.0.100\Experiment\Retina\YiKo\video workspace Br50\workspace_OU UD G5.7 5min Q85')
load(['\\192.168.0.100\Experiment\Retina\Chou\20180720 ',Dirti, ' G0',num2str(G),' 7min Br25 Q85']) %load workspace
clearvars -except  newYarray newXarray G_HMM a_data Spikes Infos name new_x new_y G Dirti

% name=['0426 HMM  ',Dirti,' G0',num2str(G),' 7min'];
% name=['0412 OU  ',Dirti,' G',num2str(G),' 5min'];
name=['test0723'];
% name=['0423 PTB HMM  ',Dirti,' G0',num2str(G),' 5min Br50'];

%% load (sorted) spikes file
% load('D:\Yiko\Files for Thesis\04122018\sortch 0412OU_UD_G19_5min_offch8')
%  load(['\\192.168.0.100\Experiment\Retina\YiKo\Experiment\04142018\New folder2\0414pos2_HMM_RL_G0',num2str(G)])  %load mat file
load('D:\Yiko\Experiment\04232018\04232018 mat\0423pos2_HMM_RL_G03_7min.mat')

%%
lumin=[];
lumin=a_data(3,:);   %Careful: can't subtract a value to the lumin series, or the correspondent Spike time would be incorrect!
figure;plot(lumin);
 
Samplingrate=20000; %fps of diode in A3
idealStimuli=[]; 
% idealStimuli=new_x;xls
idealStimuli=newXarray;


%setting for parameters
start_lum=3.44*10^4; %81*10^4;
platau_n=80;  %least number of point for platau
thre01=3.44*10^4; %cut from middle to high state
thre02=3.475*10^4; %for higher state

%if has brief pump before the video start: set a value for that
%  lumin(1:150000)=3.436*10^4;
%   lumin(6300000:end)=3.47*10^4;
%  figure;plot(lumin);

%find stimulation period
diode_start=find(lumin>=start_lum,1);
tempp=find(lumin<=start_lum);
diode_end=tempp(end); 

% plot
figure;plot(lumin)
hold on; plot(diode_start,lumin(diode_start),'r*');
hold on;plot(diode_end,lumin(diode_end),'r*');

% figure;plot(lumin(5*Samplingrate:7*Samplingrate)); hold on; plot(diode_start-5*Samplingrate,lumin(diode_start),'r*')
% figure;plot(lumin(length(lumin)-25*Samplingrate:length(lumin)-23.7*Samplingrate)); hold on; plot(diode_end-(length(lumin)-25*Samplingrate),lumin(diode_end),'r*')

%% Cut the spike time from the start and final diode timing
TimeStamps=zeros(1,2);
TimeStamps(1,1)=diode_start/Samplingrate;
TimeStamps(1,2)=diode_end/Samplingrate;
TimeStamps

 yk_spikes=[];
 for j = 1:length(Spikes)    %running through each channel
        ss = Spikes{j};
                ss(ss<TimeStamps(1,1)) = [];  %delete the spikes before TimeStamps(1)
                ss(ss>TimeStamps(1,2))=[];

        for i = 1:length(ss)
            ss(i) = ss(i)-TimeStamps(1,1);
        end
        yk_spikes{j} = ss;
 end


% delete before and after of the lumin series
lumin2=[];
lumin2=lumin;
lumin2(1:diode_start-1)=[]; %clear numbers before the diode_start
lumin2(diode_end+1-(diode_start-1):end)=[]; %clear numbers after the diode_end

%time of stimulation from diode measurement
totalTime=vpa(length(lumin2)/ Samplingrate,30);
totalTime

% assign state number 1,2,3 to three luminance state
% figure;hist(lumin,10)
numseries=[];
slumin= smooth(lumin2) ;
for jj=1:length(lumin2)
    if slumin(jj) > thre02
         numseries(jj)=3;  %highest diode value: state 3
    elseif   slumin(jj) > thre01 &&  slumin(jj)  <thre02
        numseries(jj)=2; %middle diode value: state 2
    else
        numseries(jj)=1; %lowest diode value: state 1
    end
end
dd=[]; %make the first elements = 3, since when the luminance raises up, the first  element is state 3
dd=find(numseries==3,1);
numseries(1:dd)=3;


figure; plot(lumin2)
den=0.8*(max(lumin2)-min(lumin2))/(max(numseries)-min(numseries));
numseries02=numseries*den;
hold on; plot(numseries02-mean(numseries02)+mean(lumin2)-0.2)

% find delay frames' timing
frame_time=[];  changept=[]; count=1; pt_delay=[];
changept=find(diff(numseries)~=0); %find the point that changes state
for hh=1:length(changept)-1
        frame_time(count)=length(numseries(changept(hh)+1:changept(hh+1)))/Samplingrate; %the time period for each platau
        count=count+1;
end
pt_delay=changept(find(frame_time>0.018));  %if the time of a platau is too long, it meas there is a problemed frame

figure;plot(lumin2); hold on;
plot(pt_delay,lumin2(pt_delay),'g*')
 
% for i=1:length(pt_delay)
%     figure;plot([length(lumin2(1:pt_delay(i)-Samplingrate/4)):length(lumin2(1:pt_delay(i)+Samplingrate/4))] , lumin2(pt_delay(i)-Samplingrate/4 :pt_delay(i)+Samplingrate/4)); hold on;
%     plot(pt_delay(i),lumin2(pt_delay(i)),'g*')
% end

% find transition number of point - not finished
% qwq=[];
% lumin3=diff(lumin2); lumin3(length(lumin3)+1)=0;
% qwq=find(lumin2>thre01 & lumin2<thre02 & lumin3>0);
% figure;plot(lumin2(1:10*Samplingrate));
% hold on; plot(qwq(1:50000),lumin2(qwq(1:50000)),'*')

% assign transitioning points(ex: the points from state 2 to state 3) to other states
changept=[];
changept=find(diff(numseries)~=0);
for hh=1:length(changept)-1
        if length(numseries(changept(hh)+1:changept(hh+1)))<platau_n  %to determine whether it's transitioning point
            numseries(changept(hh)+1:changept(hh+1))=numseries(changept(hh+1)+1);  %assigning the new state to transitioning point
        end
end

% % plot
figure; plot(lumin2)
recal=0.8*numseries*(max(lumin2)-min(lumin2))/(max(numseries)-min(numseries));
recal02=recal-mean(recal) +mean(lumin2)-0.02;
% hold on; plot(recal)
hold on; plot(recal02)
hold on;plot(pt_delay,lumin2(pt_delay),'g*')


for i=1:length(pt_delay)
    figure;plot([length(lumin2(1:pt_delay(i)-Samplingrate/4)):length(lumin2(1:pt_delay(i)+Samplingrate/4))] , lumin2(pt_delay(i)-Samplingrate/4 :pt_delay(i)+Samplingrate/4));
    hold on; plot([length(recal02(1:pt_delay(i)-Samplingrate/4)):length(recal02(1:pt_delay(i)+Samplingrate/4))] , recal02(pt_delay(i)-Samplingrate/4 :pt_delay(i)+Samplingrate/4))
    hold on;plot(pt_delay(i),lumin2(pt_delay(i)),'g*')
end

% assign bar position to every number in numseries(1,2,3)
stimuli_pos=[]; %the bar position value corresponding to the states of diode
Fcount=1; %counting for the idealStimuli
stimuli_pos(1)=idealStimuli(1);
skippos=zeros(1,length(numseries));

for uu=2:length(numseries)
    Ptbefore=numseries(uu-1);
    Pt=numseries(uu);
    if Pt-Ptbefore == 0  %if they are the same state
        stimuli_pos(uu)=idealStimuli(Fcount);
    elseif Pt-Ptbefore == -1 || Pt-Ptbefore==2 %normal transition
        Fcount=Fcount+1;
        stimuli_pos(uu)=idealStimuli(Fcount);
    elseif Pt-Ptbefore == 1 || Pt-Ptbefore== -2 %skipped frame
        Fcount=Fcount+2;
        skippos(uu)=1;
        stimuli_pos(uu)=idealStimuli(Fcount);
    end
end
Fcount  %the length of Fcount should be the same as the length of "idealStimuli"
pt_skip=find(skippos~=0); %find how many skipped frames

% plot
% figure; plot(numseries)
% hold on; plot(stimuli_pos*0.01-3.5);
% hold on; plot(pt_skip,numseries(pt_skip),'g*')
% hold on;plot(pt_delay,numseries(pt_delay),'g*')

% check the number of state is the same as Fcount(length of idealStimuli)
test=[];
test=find(diff(numseries)~=0);
NumFrames=length(test)+length(pt_delay)+1; %number of total different frames from diode measurement  %+1 is for the last frame
NumFrames-length(pt_delay)+length(pt_skip)==Fcount %check the assignment of bar position is correct: =1:success!
% figure; plot(numseries)
% hold on; plot(test,numseries(test),'ro')



% Previous array " stimuli_pos" is with sampling rate of mcrack, now we need to make an array with  60Hz
BinningInterval_pre=vpa(totalTime/NumFrames,40); %First: determine binning interval:
BinningInterval_pre         %symbolic form is more precise  %dont use 1/60
% BinningInterval_pre=0.01666530878596023664891020110696074195429; %turn to double > calculation is faster
digitsOld = digits(100);
BinningInterval = double(BinningInterval_pre) ;
diode_BT = [BinningInterval : BinningInterval : (diode_end-diode_start+1)/Samplingrate]; %binning time for diode
%[16,33,49,-----] %lack one frame, due to DataTime is not 150s
length(diode_BT)

bin_pos=[]; %the result we want
tog3=[];
n_1frame=BinningInterval*Samplingrate; %number of diode points in one frame
same_len_pos=[];
ratio3=[]; count=1; %ratio3 calculate the ratio for 3 states in one bin

for hj= 1:length(diode_BT) %length is same as BinningTime for Spikes %8534:8539 %8457:8461
    numT=diode_BT(hj)*Samplingrate; %the corresponding number of point to timing in diode_BT
    C=[];
    C = unique(stimuli_pos(round(numT-n_1frame)+1:round(numT)),'stable'); %unique number in one frame
%     hj
%     C
    if length(C)==1
            bin_pos(hj)=C;
    elseif length(C)==2
            pt1n=length(find(stimuli_pos(round(numT-n_1frame)+1:round(numT))==C(1))); %number of first number
            pt2n=length(find(stimuli_pos(round(numT-n_1frame)+1:round(numT))==C(2)));
                if pt1n==pt2n
                bin_pos(hj)=C(2);
                same_len_pos(hj)=round(numT-n_1frame);
                elseif pt1n > pt2n
                bin_pos(hj)=C(1);
                elseif pt1n < pt2n
                bin_pos(hj)=C(2);
               end
    elseif length(C) ==3
%             tog3(round(numT-n_1frame)+1:round(numT))=1;
            tog3(round(numT-n_1frame))=1;
            tog3(round(numT))=1;
            pt1n=length(find(stimuli_pos(round(numT-n_1frame)+1:round(numT))==C(1))); %number of first number
            pt2n=length(find(stimuli_pos(round(numT-n_1frame)+1:round(numT))==C(2)));
            pt3n=length(find(stimuli_pos(round(numT-n_1frame)+1:round(numT))==C(3)));
            err=max([pt1n, pt2n, pt3n]);
            if length(err)~=1
                hj
                break
            end
                    if err==pt1n  %assign the value that has the most numbers
                    bin_pos(hj)=C(1);
                    elseif err==pt2n
                    bin_pos(hj)=C(2);
                    elseif err==pt3n
                    bin_pos(hj)=C(3);
                    end
                    
        sumtemp=pt1n+pt2n+pt3n;
        ratio3(count,1:3)=[pt1n pt2n pt3n]./sumtemp;  %ration between these 3 values
        count=count+1;
                    
    else
         bin_pos(hj)=-5;
    end
end
length(find(tog3~=0)/2) %this is the number of points that has 3 states in one bin

disp('ratio3/ problem points: ')
ratio3(find(ratio3<0.9 & ratio3>0.1))  %3 states / if not empty, need to recalculate
length(find(same_len_pos~=0))  % 2 states / if not zero, need to recalculate


%% Saving
clearvars -except bin_pos diode_BT BinningInterval a_data Spikes yk_spikes TimeStamps  start_lum thre01 thre02 Samplingrate idealStimuli platau_n name
save(['merge_sortch_',name])




%% some plotting
% plot regions that ahs 3 states in one bin
% figure; plot(lumin2)
% recal=0.8*numseries*(max(lumin2)-min(lumin2))/(max(numseries)-min(numseries));
% recal02=recal-mean(recal) +mean(lumin2)-0.02+70;
% hold on; plot(recal02)
% % hold on;plot(pt_delay,lumin2(pt_delay),'g*')
% oao=3.49*10^4*ones(1,length(find(tog3~=0)));
% hold on; plot(find(tog3~=0),oao,'b*')

%plot
% figure; plot(lumin2)
% hold on; plot(recal02)
% hold on;plot(pt_delay,lumin2(pt_delay),'g*')
% temp=[];
% temp=mean(recal02)*ones(1,length(diode_BT));
% hold on;plot(diode_BT,temp,'b*' )
% 
% xtemp=find(tog3~=0);
% ytemp=[];
% ytemp=mean(lumin2)*ones(1,length(xtemp));
% hold on; plot(xtemp,ytemp,'ro')

% plot
% figure; autocorr(bin_pos,500)
% saveas(gcf, ['Autocorrelation bin_pos',name],'fig');
% figure; autocorr(idealStimuli,500)

% figure; plot(bin_pos);
% hold on; plot(idealStimuli));

%% transfer to ccd frame
% ccd_fps=120;
% diode_select=154431; %value before cut before and after
% 
% diodeT=(diode_select-diode_start)/Samplingrate;
% through_frame=diodeT*ccd_fps;
% ccd_startframe=249;
% targer_frame=ccd_startframe+through_frame;
% 
% targer_frame
