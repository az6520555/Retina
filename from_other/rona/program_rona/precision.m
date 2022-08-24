% Calculate precision and reliability, Rona
clear all
% close all
clc

cd('E:\google_rona\20170524\repeatHMM') ; % the folder of the files
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file) ; 
Rel = [];Prec=[];
roi = [3];
% roi = datasample(roi,1,'Replace',false);  %roi(kkk);%

for m = 1
    Nums = [];
    clearvars -except all_file n_file m roi Rel Nums Prec
    file = all_file(m).name ;
    [pathstr, name, ext] = fileparts(file);
%     directory = [pathstr,'\'];
    filename = [name,ext]
    load(filename)
    name(name=='_')='-';
    %%  TimeStamps  %%%%%%%%%%%%%%%%%%%
    if size(a_data,1)==1              %only find one analog channel, possibly cause by the setting in MC_rack
        a_data2 = a_data(1,:);
    else
        a_data2 = a_data(3,:);   
    end
    [~,locs]=findpeaks(diff(a_data2),'MINPEAKHEIGHT',20*std(diff(a_data2)));
    analog_loc = (locs)/20000;
    TimeStamps = analog_loc;

    %% spike process
    TrialSpike = [];
    starts = [TimeStamps];%
    DataTime = (TimeStamps(2)-TimeStamps(1))/2;
    BinningInterval = 2e-3; 
    BinningTime = [ 0 : BinningInterval : DataTime];

    for sweepindex = 1:length(starts)
        TimeStampsSweep = [starts(sweepindex),starts(sweepindex)+DataTime]; % forcus on which trails 
        cut_spikes = seperate_trials(Spikes,TimeStampsSweep); 
        for i = 1:length(roi)
            ii = roi(i);
            [n,xout] = hist(cut_spikes{ii},BinningTime);
            BinningSpike = n ;
            BinningSpike(:,1) = 0;BinningSpike(:,end) = 0;
            TrialSpike(sweepindex,:) = [BinningSpike];
        end      
    end
    
    %%  Plot Different Trials %%
    sweepend=sweepindex;
    figure(2);
    set(gcf,'position',[150,30,1024,900])
    h = subplot(sweepend+1,1,1);
    
    for sweepindex=1:sweepend
        subplot(sweepend+1,1,sweepindex);
%         if size(TrialSpike)>1
%             plot(BinningTime,TrialSpike{1}(sweepindex,:),'b',BinningTime,TrialSpike{2}(sweepindex,:),'g');
%         else
            plot(BinningTime,TrialSpike(sweepindex,:),'b');
%         end
    end     
    subplot(sweepend+1,1,sweepindex+1);
    y1 = a_data(3,TimeStamps(1)*20000:(TimeStamps(1)+DataTime)*20000);
    x1 = 1/20000:1/20000:length(y1)/20000;
    plot(x1,y1,'r');
    samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
    set(get(h,'title'),'string',[name,'  ch',num2str(roi)]);
    ylim([4.6*10^4,4.9*10^4]);
    TrialSpike(find(TrialSpike~=0))=1;
    
    figure;imagesc(BinningTime,[1:sweepend],TrialSpike);
    title([name,'  ch',num2str(roi)]);  
    xlabel('t(s)');ylabel('trial');
    colormap gray;
    cmap = colormap(gray); cmap = cmap(end:-1:1, :); colormap(cmap);
        %% precision & reliability
    window = 5; %tolerence=window*BinningInterval
    rel = []; 
    prec = []; 
    for i = 1:window:size(TrialSpike,2)-(window-1)  %discrete; moving window:i = 1:size(TrialSpike,2)-(window-1)
        ws =  TrialSpike(:,i:i+window-1);  % window spike
        sws = sum(TrialSpike(:,i:i+window-1),2); % sum of window spike
        sws(find(sws>=1)) = 1;
        [a,b] = find(ws~=0);
        if isempty(b)~=1
            prec = [prec std(b)*BinningInterval*10^3];
            rel = [rel; sum(sws)/size(TrialSpike,1)];
        end
    end
    Prec = [Prec; mean(prec)];
    Rel = [Rel; mean(rel)];
    
%     %% precision
%     window = 2; %tolerence=window*BinningInterval
%     rel = []; 
%     prec = []; 
%     for i = 1:size(TrialSpike,2)-(window-1)  %round(size(KKK,2)/window)-window  %
%     %     temp = sum(TrialSpike(:,i:i+window),2);
%         temp =  TrialSpike(:,i:i+window-1);  %sum(KKK(:,window*i:(i+1)*window),1);
%         [a,b] = find(temp~=0);
%         if isempty(b)~=1
%             rel = [rel;length(b)/length(temp)];
%     %         nums = [nums; temp(b)];
%             prec = [prec std(b)*BinningInterval*10^3];
%         end
%     end
%     Prec = [Prec; mean(prec)];
% % figure; hist(prec)
%     %% reliability
%     for i = 1:size(TrialSpike,2)-(window-1)  %round(size(KKK,2)/window)-window  %
%         temp = sum(TrialSpike(:,i:i+window-1),1);
%         [a,b] = find(temp~=0);
%         if isempty(b)~=1
%             rel = [rel; length(b)/length(temp)];
%         end
%     end
%     Rel = [Rel; mean(rel)];
% 
%     %% firing rate
%     fr=mean(mean(TrialSpike,2))/BinningInterval;
end
% plot(1./([3 5 10 30]),Prec,'-o')
% figure;plotyy(1./([3 5 10 30]),Rel,1./([3 5 10 30]),Prec)
% plot(32./([5 10 20 40 80 200]),Prec,'-o')
% plot(32./([5 10 20 40 80 200]),Rel,'-o')
% hold on



