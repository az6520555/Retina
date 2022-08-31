clear
close all
Tot = 300;  
dt = 0.01; 
T = dt:dt:Tot;
% tau=0.5;
D =4; %dynamical range
at = 60;
a2=[];ey=[];ey1=[];
saving=0;
tau2=[0.5 1]
c={'k','r'}
for jj=1:2
    tau=tau2(jj)
    L = zeros(1,length(T));
    for i = 1:length(T)-1
        L(i+1) = L(i) + (-L(i)/tau + randn*sqrt(D/dt))*dt;
    end


%     path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\stimulus_data\OU_different_contrast\20210716\';
%     figure(1);plot(T,L,'linewidth',1);xlabel('time (s)')
    cutoff=1;
    corrtimef=[];

    [b,a]=butter(2,2*cutoff*dt,'low')
    Lf=filter(b,a,L);

    [val,lags] =autocorr(L,300);
    [M,I]=max(val);
    % figure(1);plot(lags,val)
    Lfshift=Lf(-lags(I)+1:end);

    figure(1);hold on
    plot(lags,val,'linewidth',1,'color',c{jj})
end
box on
ylim([0,1])
title('Autocorrelation of OU')
xlabel('number of lags')
legend(('\tau = 0.5 s'),('\tau = 1 s'))
set(gcf,'Position',[300,300,300,300])
% figure(1);hold on
% plot(T(1:end+lags(I)),Lfshift,'linewidth',1)

% legend(('OU'),('LPOU f_c=1Hz'))
%% generate stimulus
% figure(4);hold on;
% [t,ey,a2]=generate_stimulus(L,10,2,dt);  % (signal, mean intensity, times of std)
% plot(t,ey,'linewidth',1);
% if saving
%     save([path,'OU_tau=',num2str(tau*1000),'ms_cutoff=0_mean10_amp2.mat'],'ey','a2','t')
% end
% [t,ey,a2]=generate_stimulus(Lfshift,10,2,dt);  % (signal, mean intensity, times of std)
% plot(t,ey,'linewidth',1);
% if saving
%     save([path,'OU_tau=',num2str(tau*1000),'ms_cutoff=1_mean10_amp2.mat'],'ey','a2','t')
% end
% 
% % different mean same contrast
% figure(78);hold on
% inten_coeff=[1.3,0.6,0.2];
% for i=1:3
%     [t,ey,a2]=generate_stimulus(Lfshift,10,2,dt);
%     ey=inten_coeff(i)*ey;
% %     plot(t,ey)
%     if saving
%        save([path,'OU_tau=',num2str(tau*1000),'ms_cutoff=1_coeff=',num2str(inten_coeff(i)),'.mat'],'ey','a2','t')
%     end
% end
% 
% % different mean same amplitude
% figure(87);hold on
% % mean_inten_set=[15,10,5];
% amp_set=[0.5,1,1.5,3];
% for i=1:length(amp_set)
%     mean=10;
%     amp=amp_set(i);
%     [t,ey,a2]=generate_stimulus(Lfshift,mean,amp,dt);
% %     plot(t,ey)
%     if saving
%       save([path,'OU_tau=',num2str(tau*1000),'ms_cutoff=1_mean=',num2str(mean),'_amp=',num2str(amp),'.mat'],'ey','a2','t')
%     end
% end
% 
% 
% % ======== thesis setting ======
% figure(4)
% box on
% title('OU and LPOU')
% xlabel('time (s)')
% ylabel('light intensity (mW/m^2)')
% set(gcf,'Position',[300,300,600,300])
% xlim([100,104])
% legend(('OU'),('LPOU with f_c=1 Hz'))