clear
close all
cd('F:\§Úªº¶³ºİµwºĞ\Retina exp\exp data\Spatial stimuli\20210504\data')
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
rate=20000;
bin=0.001;
files=[45];

for z = 1:length(files)
    file = all_file(files(z)).name;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    stimulus=a_data(3,:)-min(a_data(3,:));
    stimulus=stimulus/max(stimulus);
    t=1/rate:1/rate:length(a_data(3,:))/rate;
    BinningTime=bin:bin:size(a_data,2)/rate;
    
%     for i = 1:60
%         [n,time_bin]=hist(Spikes{i},BinningTime);
%         w = gausswin(100);
%         n_fil=filter(w,1,n);
%         if length(Spikes{i})>100
%             figure(i);hold on
%             yyaxis left
%             plot(t,stimulus)
%             yyaxis right
%             plot(time_bin,n_fil)
%             plot(time_bin,n,'color','g')
%         end
%     end
    % ======= plot peak position of gaussian filtered spike ========
    sti_sm=smooth(stimulus,100);
    sti_diff=diff(stimulus);
    threshold=0.5;
    sti_near_zero=-abs(stimulus-threshold);
    [~,half_locs]=findpeaks(sti_near_zero,'MINPEAKDISTANCE',100,'MINPEAKHEIGHT',-0.2);
    pulse_start_locs=half_locs(1:2:end);
    pulse_end_locs=half_locs(2:2:end);
    figure(87); set(gcf,'Name','stimulus');
    plot(t,stimulus);hold on
    plot(t(pulse_start_locs),0.5*ones(1,length(pulse_start_locs)),'o')
    plot(t(pulse_end_locs),0.5*ones(1,length(pulse_end_locs)),'go')
    
    % ====== Sum the spikes in a single flicker (using hist) ==========
    total_pulse=length(pulse_start_locs);
    adap_time=10; % eliminate the duration of adaptation
    n_select_pulses=find(t(pulse_start_locs)>t(pulse_start_locs(1))+adap_time);
    for channel=1:60
        if length(Spikes{channel})>300
            [n_spikes n_pulse]=hist(Spikes{channel},t(pulse_start_locs(n_select_pulses)));
            figure;
            set(gcf,'Position',[10 50 1000 950])
            subplot(2,1,1)
            plot(n_pulse(2:end-1),n_spikes(2:end-1),'-o')
            ylim([0 max(n_spikes(2:end-1))+1])
            subplot(2,1,2)
            plot(n_spikes(2:end-1),n_spikes(3:end),'.')
            title(['channel ',num2str(channel),',  spike number=',num2str(length(Spikes{channel}))])
        end
    end
    
    
end
