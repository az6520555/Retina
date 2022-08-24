%% Gaussian pulse
dt=0.1;
% tau=0.2;
% xp=-tau*3/dt:tau*3/dt;
% ey1=[];
% ey1=exp(-(xp*dt/tau).^2);
x=-100:dt:100;
ey1=sin(x);
%%
m=20;
b=(3+m)/2;
x=ey1;
FBsum=zeros(1,length(x));
for k=0:m-1
    fb=[];
    for i=1:length(x)
        ck=(k+1)/m;
        ytime=i/dt-(m-k)*dt;
        ind=fix(ytime*dt);
        if ind>0
            fb(i)=ck*x(ind);
        elseif ind<=0
            fb(i)=0;
        end
    end
    FBsum=FBsum+fb;
end
y=(3+m)/2*x-FBsum;
t=1*dt:1*dt:length(x)*dt;
figure(1);plot(t,x,t,y,'linewidth',2)
xlim([0 max(t)])
xlabel('time (s)')
legend('Input','Output (m=150)')