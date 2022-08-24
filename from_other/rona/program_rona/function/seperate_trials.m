function[spikes] = seperate_trials(spikes,starts,clear_before_first_trial)
    if ~exist('clear_before_first_trial')
        clear_before_first_trial = 1;
    end
    for j = 1:length(spikes)         %running through each channel
        ss = spikes{j}; %remove spikes before and after the stimulus
            if clear_before_first_trial == 1
                ss(ss<starts(1)) = [];
                ss(ss>starts(end))=[];
                ss(ss==0)=[];
            end
        for i = 1:length(ss) % transfer the spike times to the form related to the beginning of the trial
            loc = find(ss(i)>starts);
            if isempty(loc);  continue;   end;
            ss(i) = ss(i)-starts(loc(end));
        end
    %     ss = ss -1;
        spikes{j} = ss;
    end
end