%compare the spike trial of same stimuli
clear
close all
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210504\data')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
k=0;
rest=2;
np=50; %number of repeat of same stimuli period
loop_ind=[45];
fig_num=0;
rate=20000;
fig=zeros(1,60);

for z = loop_ind
    
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA roi fignames date path_sort k ...
        k rest1 rest np fig_num rate ax1 ax2 fig

    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    % =========deal with redundant stimuli==========
    Indt=1;
    for i=1:length(TimeStamps)-1
        Dt=TimeStamps(i+1)-TimeStamps(i);
        if Dt>rest-1
            Indt=Indt+1;
            redun=TimeStamps(i+1);
        end
    end
    if Indt>np
        TimeStamps=TimeStamps(TimeStamps<redun);
    end
    
    TimeStamps=TimeStamps(2:end); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ======================================
    
    sti=a_data(3,:);
    sti_t=0:1/rate:1/rate*(length(sti)-1);
    indtt=0;
    tt_st(1)=TimeStamps(1);
    for i=1:length(TimeStamps)-1
        tt=TimeStamps(i+1)-TimeStamps(i);
        if tt>rest-1
            indtt=indtt+1;
            tt_st(indtt+1)=TimeStamps(i+1);
            tt_end(indtt)=TimeStamps(i);
        end
    end
    tt_end(end+1)=TimeStamps(end);

    ind_sti=1;
    [Spikes_re,TimeStamps2]=seperate_trials2(Spikes,TimeStamps,rest,np); % seperate trials of different period 
    for a = 1:size(Spikes_re,1)
        fig_num=fig_num+1;
        sp=Spikes_re(a,:);

        k=0;
        
        y=sti(sti_t>tt_st(ind_sti) & sti_t<tt_end(ind_sti));
        x=sti_t(sti_t>tt_st(ind_sti) & sti_t<tt_end(ind_sti))-tt_st(ind_sti);
            
        ind_sti=ind_sti+1; % index of stimuli

        k=k+1; % subplot number

        % === combine all channels and plot firing rate ===
%         sp1_all=cat(2,sp{:});
        BinningInterval=0.010;
        BinningTime=0:BinningInterval:TimeStamps2{a}(end);
%         [n,SpikeTime]=hist(sp1_all,BinningTime);
%         FR=n/BinningInterval;     
%         sn=smooth(n);

%         figure(fig_num*10)
%         subplot(np+1,1,k)
%         plot(SpikeTime,n);hold on
%         plot(SpikeTime,sn);hold on % plot smoothed firing rate

%         mpd=(TimeStamps2{a}(2)-TimeStamps2{a}(1))*(1/BinningInterval)/2;
%         [pks,locs]=findpeaks(sn,'MinPeakHeight',1,'MinPeakDistance',mpd); %find peaks of smoothed firing rate
%         plot(SpikeTime(locs),pks,'diamond')
%         ylim([0 200]);
%             
%         figure(fig_num*10)
%         subplot(np+1,1,np+1)
%         plot(x,y,'r')
%         samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
    end
    
    % =========== single cell summation of the trials ===============
    for channel = 1:60 % selected cell
        
        spikes_cell=Spikes_re(:,channel);
        trial_sum = cat(2,spikes_cell{:});
        [n_cell,SpikeTime_cell]=hist(trial_sum,BinningTime);
        FR_cell=n_cell/np/BinningInterval;
        sFR=smooth(FR_cell);
        if fig(channel) == 0
            fig(channel)=figure(channel);title(name)
            ax1(channel)=subplot(2,1,1);hold on;
            ax2(channel)=subplot(2,1,2);hold on
        end
        
        plot(ax1(channel),SpikeTime_cell,FR_cell);
        lgd=legend(ax1(channel),['9hz',' channel',num2str(channel)])
        lgd.FontSize = 50;
  
%         plot(SpikeTime_cell,sFR)
        plot(ax2(channel),x,y,'LineWidth',2);
        samexaxis('abc','xmt','on','ytac','join','yld',1);
    end
end

