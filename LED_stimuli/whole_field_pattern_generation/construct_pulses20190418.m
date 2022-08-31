%% cosntruct stimuli (fixed pulse width)
rate = 20000;
f_start = 3;
f_end = 18;
p_start = 1/f_start;
p_end = 1/f_end;
pw = 0.05; %pulse width
t_trial = 30; % time of a single trial
n_trial = 5; % number of trials of same period
p = linspace(p_start,p_end,20);
rest = t_trial;
rest2 = 150;

ey = [];
for i = length(p)/2+1:length(p)
    np = t_trial/p(i); % number of pulse in a trial
    A = [];
    B = [];
    singlepulse = [ones(1,round(pw*rate)) zeros(1,round((p(i)-pw)*rate))];
    for j = 1:np
        A = [A singlepulse];
    end
    for k = 1:n_trial
        B = [B A zeros(1,round(rest*rate))];
    end
    ey = [ey B zeros(1,round(rest2*rate))];
end

ey = 0.15*ey;
a2 = ey;
t=1/rate:1/rate:length(ey)/rate;
plot(t,ey)

save('F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\stimulus data\constant_period_pulse_20190418\20190418_pulses_2.mat','ey','a2','t')