%compare the spike trial of same stimuli
clear
close all
cd('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\temporal pattern\20190419')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
fignames={all_file.name}; % assign name for figures
k=0;
rest1=30;
rest2=150;
np1=5; %number of repeat of same stimuli period
np2=10; % number of different stimuli period
loop_ind=[15 14];
fig_num=0;
rate=20000;

for z = loop_ind
    
    clearvars -except all_file n_file z SamplingRate cc ey isi2 statispks statistime w fr information rr STA roi fignames date path_sort k ...
        k rest1 rest2 np1 np2 fig_num rate

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
        if Dt>rest2-1
            Indt=Indt+1;
            redun=TimeStamps(i+1);
        end
    end
    if Indt>np2
        TimeStamps=TimeStamps(TimeStamps<redun);
    end
    % ======================================
    
    sti=a_data(3,:);
    sti_t=0:1/rate:1/rate*(length(sti)-1);
    indtt=0;
    tt_st(1)=TimeStamps(1);
    for i=1:length(TimeStamps)-1
        tt=TimeStamps(i+1)-TimeStamps(i);
        if tt>rest1-1
            indtt=indtt+1;
            tt_st(indtt+1)=TimeStamps(i+1);
            tt_end(indtt)=TimeStamps(i);
        end
    end
    tt_end(end+1)=TimeStamps(end);

    ind_sti=1;
    [Spikes_re,TimeStamps2]=seperate_trials2(Spikes,TimeStamps,rest2,np2); % seperate trials of different period 
    for a = 1:size(Spikes_re,1)
        fig_num=fig_num+1;
        sp=Spikes_re(a,:);
        [Spikes_sing,TimeStamps3]=seperate_trials2(sp,TimeStamps2{a},rest1,np1); % seperate repeat trials
        
        k=0;
        
        y=sti(sti_t>tt_st(ind_sti) & sti_t<tt_end(ind_sti));
        x=sti_t(sti_t>tt_st(ind_sti) & sti_t<tt_end(ind_sti))-tt_st(ind_sti);
            
        for j=1:size(Spikes_sing,1)
            
            
            ind_sti=ind_sti+1; % index of stimuli
            sp1=Spikes_sing(j,:);
            
            k=k+1; % subplot number
            
            % === combine all channels and plot firing rate ===
            sp1_all=cat(2,sp1{:});
            BinningInterval=0.01;
            BinningTime=[0:BinningInterval:TimeStamps3{j}(end)];
            [n,SpikeTime]=hist(sp1_all,BinningTime);
%             FR=n/BinningInterval;     
            sn=smooth(n);

            figure(fig_num*10)
            subplot(np1+1,1,k)
            plot(SpikeTime,n);hold on
            plot(SpikeTime,sn);hold on % plot smoothed firing rate
            
            
            mpd=(TimeStamps3{j}(2)-TimeStamps3{j}(1))*(1/BinningInterval)/2;
            [pks,locs]=findpeaks(sn,'MinPeakHeight',1,'MinPeakDistance',mpd); %find peaks of smoothed firing rate
            plot(SpikeTime(locs),pks,'diamond')
            ylim([0 15]);
            
            % === raster plot ===
%             for i=1:60
%                 if isempty(sp1{i})==1
%                     sp1{i}=0;
%                 end
%             end
            
%             figure(fig_num);
%             subplot(np1+1,1,k)
%             LineFormat.Color = [0.3 0.3 0.3];
%             plotSpikeRaster(sp1,'PlotType','vertline','LineFormat',LineFormat);
%             hold on
        end
%         figure(fig_num);
%         subplot(np1+1,1,np1+1)
%         plot(x,y,'r')
%         samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
        
        figure(fig_num*10)
        subplot(np1+1,1,np1+1)
        plot(x,y,'r')
        samexaxis('abc','xmt','on','ytac','join','yld',1) % function need ''parseArgs'' function, from matlab center
    end
end


    
