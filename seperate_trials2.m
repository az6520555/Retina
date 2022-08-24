% ========= use after removing the redundant stimuli =========
function [NewSpikes,NewTimeStamps] = seperate_trials2(SPIKES,TIMESTAMPS,rest,NP)
    indt=0; % the nth trial  
    trial_start(1)=TIMESTAMPS(1);
    ind_start(1)=1;
    for m=1:length(TIMESTAMPS)-1
        dt=TIMESTAMPS(m+1)-TIMESTAMPS(m);
        if dt>rest-1
            indt=indt+1;
            ind_start(indt+1)=m+1;
            ind_end(indt)=m;
            trial_start(indt+1)=TIMESTAMPS(m+1);
            trial_end(indt)=TIMESTAMPS(m);
        end
    end
    ind_end(end+1)=length(TIMESTAMPS);
    trial_end(end+1)=TIMESTAMPS(end);
    
    for m=1:length(SPIKES)
        if length(trial_start)>NP
            disp('deal with the TimeStamps, OK?');
            return
        else
            for n=1:length(trial_start)
                ss=[];
                ss=SPIKES{m};
                ss=ss(ss>trial_start(n) & ss<trial_end(n))-trial_start(n);
                NewSpikes{n,m}=ss;
            end
        end
    end
    
    for p=1:length(trial_start)
        NewTimeStamps{p}=TIMESTAMPS(ind_start(p):ind_end(p))-TIMESTAMPS(ind_start(p));
    end
end
    