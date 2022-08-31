% test stimuli
% LED_pattern(stimu,filename,pw)

clear
close all

parameters=...
    {'cn',0,0;...
    'oo',0,0;...
    'gn',0,0.3;...
    'gn',0,0.1;...
    'gn',0,0.6;...
    'gf',0,0.6;...
    'gf',0,0.3;...
    'gf',0,0.1;...
    'jk',0,0};
    
    
stimu_set=parameters(:,1);  %
filename_set=parameters(:,2);
pw_set=parameters(:,3);

ey_total=[];
a2_total=[];
for i=1:length(stimu_set)
    stimu=stimu_set{i};
    filename=filename_set{i};
    pw=pw_set{i};
    
    [ey,a2]=LED_pattern(stimu,filename,pw);
    
    ey_total=[ey_total ey zeros(1,300*20000)];
    a2_total=[a2_total a2 zeros(1,300*20000)];
    
end
figure;subplot(2,1,1)
plot(ey_total);
subplot(2,1,2)
plot(a2_total)

% for i=1:length(stimu_set)
%     stimu=stimu_set{i};
%     filename=filename_set{i};
%     pw=pw_set{i};
%     
%     DAQ_LED_autotest2
%     
%     pause(300)
% end
