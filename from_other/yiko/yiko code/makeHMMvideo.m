%% calculate spatial error
% load workspace from Oct_calibration  :  load file: workspace_spatial_cal.mat
% clearvars -except  dotPositionMatrix detect_pt ideal_pt
% 
% mea_size=503;
% mea_size_bm=541; %bigger mea size , from luminance calibrated region
% meaCenter_x=630; 
% meaCenter_y=573; 
% leftx_bm=meaCenter_x-(mea_size_bm-1)/2; %the first x position of the bigger mea region(luminance calibrated region) on LED screen
% lefty_bm=meaCenter_y-(mea_size_bm-1)/2;
% bar_le=(mea_size-1)/2; %half of bar length / pixel number on LCD /total length = mea_size = 1919 um
% bar_wid=11; 
% 
% x_error=[];
% y_error=[];
% for i=1:8
%     for j=1:8
%         joko=cell2mat(detect_pt(i,j))-cell2mat(ideal_pt(i,j));
%         x_error(i,j)=joko(1);
%         y_error(i,j)=joko(2);
%     end
% end
% 
% %since ideal and detect_pt are point of ccd pixel > first transform to um
% pixelccd=2.24; %the size(um) of 1 pixel of ccd
% x_error=x_error.*pixelccd; %turn it to um
% y_error=y_error.*pixelccd;
% 
% %then use
% pixelscr=(2.87+2.816)/2; %the size(um) of 1 pixel of screen
% x_error(abs(x_error)<pixelscr)=0;
% y_error(abs(y_error)<pixelscr)=0;
% 
% x_error02=floor(x_error./pixelscr); %how many pixel on screen need to be adjusted
% y_error02=floor(y_error./pixelscr);
% 
% %minus one pixel of screen
% x_error03=x_error02;
% y_error03=y_error02;
% for k=1:64
%     if x_error02(k)<0
%         x_error03(k)=x_error02(k)+1;
%     elseif x_error02(k)>0
%         x_error03(k)=x_error02(k)-1;
%     end
%     if y_error02(k)<0
%         y_error03(k)=y_error02(k)+1;
%     elseif y_error02(k)>0
%         y_error03(k)=y_error02(k)-1;
%     end
%     
% end
% 
% 
% %
% y_error03
% 0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% -1	0	0	0	0	0	0	0
% -2	-1	-1	-1	0	-1	0	-1
% 
% x_error03
% 1	0	0	0	0	0	0	0
% 1	0	0	0	0	0	0	0
% 1	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% 0	0	0	0	0	0	0	0
% 1	0	0	0	0	0	0	0
% 1	0	0	0	0	0	0	0



%% HMM RL
clear all
G_list=[ 2.5 3 4.3 4.5 5.3 6.3 6.5 7.5 9 12 20];  %list of Gamma valau
countt=1;
for Gvalue=G_list
clearvars -except G_list Gvalue countt
%load the workspace of the stimualtion, containing: channel_map, clean_LED and dotPositionMatrix
load('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\workspace_for_video') 


%for randon number files ( I specifically choose some certain random seed series
cd('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\rn_workspace') ;  
all_file = dir('*.mat') ;
file = all_file(countt).name ;
[pathstr, name, ext] = fileparts(file);
directory = [pathstr,'\'];
filename = [name,ext];
load([filename]);
name=[name];
name
Gvalue
countt=countt+1;

fps =60;  %freq of the screen flipping 
T=7*60; %second
dt=1/fps;
T=dt:dt:T;

G_HMM =Gvalue; % damping / only G will influence correlation time
D_HMM = 2700000; %dynamical range
omega =G_HMM/2.12;   % omega = G/(2w)=1.06; follow Bielak's overdamped dynamics/ 2015PNAS

%R-L
init=channel_map{68}; %L
final=channel_map{69}; %R

Xarray = zeros(1,length(T));
Xarray(1,1)=0; % since the mean value of damped eq is zero
Vx = zeros(1,length(T));
%Use rntest(t)!!!
for t = 1:length(T)-1
        Xarray(t+1) = Xarray(t) + Vx(t)*dt;
        Vx(t+1) = (1-G_HMM*dt)*Vx(t) - omega^2*Xarray(t)*dt + sqrt(dt*D_HMM)*rntest(t); 
end
% Normalize to proper moving range
nrx=abs((final(1)-init(1))/(max(Xarray)-min(Xarray)));
Xarray2=Xarray*nrx;
mdist=abs(min(Xarray2)-init(1)); %rearrange the boundary values
    if min(Xarray2)>init(1)
     mdist=-mdist;
    end
Xarray3=Xarray2+mdist;
newXarray=round(Xarray3); 
Y =round((init(2)+final(2))/2); 

% clean_LED=1.2*clean_LED./255; %before 04022018
clean_LED=clean_LED./255; %make it to 0-1 range for double (due to imwrite format)
clean_LED(clean_LED>1)=1;

mea_size=503;
mea_size_bm=541; %bigger mea size , from luminance calibrated region
meaCenter_x=630; 
meaCenter_y=573; 
leftx_bm=meaCenter_x-(mea_size_bm-1)/2; %the first x position of the bigger mea region(luminance calibrated region) on LED screen
lefty_bm=meaCenter_y-(mea_size_bm-1)/2;

bar_le=(mea_size-1)/2; %half of bar length / pixel number on LCD /total length = mea_size = 1919 um
bar_wid=11; %half of bar width / total length = 11*2+1=23 pixels = 65 um

%rearrange dotPositionMatrix format
dotPosition_x=[];   dotPosition_y=[];
dotPosition_x=dotPositionMatrix(1,1:8:end);
dotPosition_y=dotPositionMatrix(2,1:8);

number_repeat=repmat([1:12000],1,3);

%number file
cd('D:\Yiko\Timing stuff\numberFlip')
all_fileNumber= dir('*.jpeg') ;
%video frame file
name=['0422 HMM RL G',num2str(G_HMM) ,' 7min Br25 Q85'];
name
% ori_dir =['D:\Yiko\HMM video frame\',name];   

%video setting
Time=T; %sec
video_fps=fps;
writerObj = VideoWriter(['\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\0422 HMM video 7min Br25 Q85\', name,'.avi']);  %change video name here!
writerObj.FrameRate = video_fps;
writerObj.Quality = 85;
open(writerObj);
%start part: dark adaptation
for mm=1:60
img=0*ones(1024,1280);   
writeVideo(writerObj,  img);
end


for kk =1:length(T)
a=zeros(1024,1280);%full screen pixel matrix %it's the LED screen size

%HMM RL bar trajectory
X=newXarray(kk);

if X-bar_wid <= dotPosition_x(1) %the upper-left corner of the mea region
    ypart=(Y+bar_le-(Y-bar_le))/14;  %divide the total length into 14 parts
    portion01=round(ypart*5); %according to x_error03
    portion02=round(ypart*11);
    portion03=round(ypart*13);
%part1
    py1=Y-bar_le;
    py2=Y-bar_le+portion01;
    px1=X-bar_wid-1;
    px2= X+bar_wid;
    a(py1:py2, px1:px2)=clean_LED(py1-round(lefty_bm):py2-round(lefty_bm),  px1-round(leftx_bm):px2-round(leftx_bm)); 
%part2
    py3=Y-bar_le+portion01+1;
    py4=Y-bar_le+portion02;
    px3=X-bar_wid;
    px4= X+bar_wid;
    a(py3:py4 ,  px3:px4 )=clean_LED(py3-round(lefty_bm):py4-round(lefty_bm),  px3-round(leftx_bm):px4-round(leftx_bm)); 
%part3
    py5=Y-bar_le+portion02+1;
    py6=Y-bar_le+portion03;
    px5=X-bar_wid-1 ;
    px6= X+bar_wid;
    a(py5:py6 ,  px5:px6 )=clean_LED(py5-round(lefty_bm):py6-round(lefty_bm),  px5-round(leftx_bm):px6-round(leftx_bm)); 
%part4
    py7=Y-bar_le+portion03+1;
    py8=Y+bar_le;
    px7=X-bar_wid-1;
    px8= X+bar_wid;
    ptmiddle=round((dotPosition_x(2)+dotPosition_x(1))/2);
    a(py7:py8 , px7:px8 )=clean_LED(py7-round(lefty_bm):py8-round(lefty_bm),  px7-round(leftx_bm):px8-round(leftx_bm));
        if px7 < dotPosition_x(1)
                 a(py8:py8+2 ,  px7:dotPosition_x(1))=clean_LED(py8-round(lefty_bm):py8+2-round(lefty_bm),  px7-round(leftx_bm): dotPosition_x(1)-round(leftx_bm));
        end
        
else
barX=X-round(leftx_bm);
barY=round(Y)-round(lefty_bm);
barImage=clean_LED(barY-bar_le: barY+bar_le, barX-bar_wid: barX+bar_wid); 
a(Y-bar_le:Y+bar_le,X-bar_wid:X+bar_wid)=barImage;
end


%number representing bar position
imgNum=imread(all_fileNumber(number_repeat(kk)).name); %read in image of number
imgNum2=rgb2gray(imgNum);
% figure;imshow(imgNum2)
%sharpen the number image edge
imgNum2(imgNum2<35)=0;
imgNum2(imgNum2>=35)=255;
imgNum2(imgNum2==255)=1;
temp1=420; %position of the number image
a(temp1+1:temp1+43, 1232: 1280)=imgNum2;

%square_flicker
if mod(kk,3)==1 %odd number
a(500-35:500+35,1230:1280)=1; % white square
elseif mod(kk,3)==2
a(500-35:500+35,1230:1280)=0.2; %gray 
else
a(500-35:500+35,1230:1280)=0; % dark
end

% imwrite(a, [ori_dir, '\', num2str(kk,'%06i'),'.png']);
%imwrite conversion would change data value and class/ Indexed images are converted to RGB before writing out JPEG files

img=[];
img=a;
writeVideo(writerObj,  img);
end

%end part video
for mm=1:10
            img=0*ones(1024,1280);
            img(500-35:500+35,1230:1280)=0.2; %gray 
            writeVideo(writerObj,  img);
end
close(writerObj);

%save parameters needed
ori_dir02 =['D:\Yiko\HMM video frame\workspace_',name,'.mat'];   
save(ori_dir02)


end



%% HMM UD
clear all
G_list=[ 2.5 3 4.3 4.5 5.3 6.3 6.5 7.5 9 12 20];
countt=1;
for Gvalue=G_list
clearvars -except G_list Gvalue countt
load('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\workspace_for_video')

%for randon number files
cd('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\rn_workspace') ;  
all_file = dir('*.mat') ;
file = all_file(countt).name ;
[pathstr, name, ext] = fileparts(file);
directory = [pathstr,'\'];
filename = [name,ext];
load([filename]);
name=[name];
name
Gvalue
countt=countt+1;


fps =60;  %freq of the screen flipping 
T=7*60; %second
dt=1/fps;
T=dt:dt:T;

G_HMM = Gvalue; % damping / only G will influence correlation time
D_HMM = 2700000; %dynamical range
omega =G_HMM/2.12;   % omega = G/(2w)=1.06; follow Bielak's overdamped dynamics/ 2015PNAS

%U-D
init=channel_map{66}; %U
final=channel_map{67}; %D

Yarray = zeros(1,length(T));
Yarray(1,1)=0; % since the mean value of damped eq is zero
Vy = zeros(1,length(T));
for t = 1:length(T)-1
        Yarray(t+1) = Yarray(t) + Vy(t)*dt;  
        Vy(t+1) = (1-G_HMM*dt)*Vy(t) - omega^2*Yarray(t)*dt + sqrt(dt*D_HMM)*rntest(t); % rntest(t): gaussian noise with zero mean and unit variance
end
%Normalize to proper moving range
nry=(final(2)-init(2))/(max(Yarray)-min(Yarray));
Yarray2=Yarray*nry;
mdist=abs(min(Yarray2)-init(2));
    if min(Yarray2)>init(2)
     mdist=-mdist;
    end
Yarray3=Yarray2+mdist;
newYarray=round(Yarray3); 
X =round((init(1)+final(1))/2); 

clean_LED=clean_LED./255;%make it to 0-1 range for double (due to imwrite format)
clean_LED(clean_LED>1)=1;

mea_size=503;
mea_size_bm=541; %bigger mea size , from luminance calibrated region
meaCenter_x=630; 
meaCenter_y=573; 
leftx_bm=meaCenter_x-(mea_size_bm-1)/2; %the first x position of the bigger mea region(luminance calibrated region) on LED screen
lefty_bm=meaCenter_y-(mea_size_bm-1)/2;

bar_le=(mea_size-1)/2; %half of bar length / pixel number on LCD /total length = mea_size = 1919 um
bar_wid=11; %half of bar width / total length = 23 pixels = 65 um

%rearrange dotPositionMatrix format
dotPosition_x=[];   dotPosition_y=[];
dotPosition_x=dotPositionMatrix(1,1:8:end);
dotPosition_y=dotPositionMatrix(2,1:8);

number_repeat=repmat([1:12000],1,3);

%number file
cd('D:\Yiko\Timing stuff\numberFlip')
all_fileNumber= dir('*.jpeg') ;
%video frame file
name=['HMM UD G',num2str(G_HMM) ,' 10min Q85'];
name
% ori_dir =['D:\Yiko\HMM video frame\',name];   

%video setting
Time=T; %sec
video_fps=fps;
writerObj = VideoWriter(['D:\Yiko\HMM video frame\', name,'.avi']);  %change video name here!
writerObj.FrameRate = video_fps;
writerObj.Quality = 85;
open(writerObj);
%start part: dark adaptation
for kk=1:60
img=0*ones(1024,1280);   
writeVideo(writerObj,  img);
end


for kk =1:length(T)
a=zeros(1024,1280);%full screen pixel matrix %it's the LED screen size
%HMM UD bar trajectory
Y=newYarray(kk);

if Y+bar_wid >= dotPosition_y(end)
    xpart=(X+bar_le-(X-bar_le))/14;
    portion01=round(xpart);
    portion02=round(xpart*7);
    portion03=round(xpart*9);
    portion04=round(xpart*11);
    portion05=round(xpart*13);
%part1
    py1=Y-bar_wid;
    py2=Y+bar_wid+2;
    px1=X-bar_le;
    px2=X-bar_le+portion01;
    a(py1:py2, px1:px2)=clean_LED(py1-round(lefty_bm):py2-round(lefty_bm),  px1-round(leftx_bm):px2-round(leftx_bm)); 
%part2
    py3=Y-bar_wid;
    py4=Y+bar_wid+1;
    px3=X-bar_le+portion01+1;
    px4=X-bar_le+portion02;
    a(py3:py4, px3:px4 )=clean_LED(py3-round(lefty_bm):py4-round(lefty_bm),  px3-round(leftx_bm):px4-round(leftx_bm)); 
%part3
    py5=Y-bar_wid;
    py6=Y+bar_wid;
    px5=X-bar_le+portion02+1;
    px6=X-bar_le+portion03;
    a(py5:py6, px5:px6 )=clean_LED(py5-round(lefty_bm):py6-round(lefty_bm),  px5-round(leftx_bm):px6-round(leftx_bm)); 
%part4
    py7=Y-bar_wid;
    py8=Y+bar_wid+1;
    px7=X-bar_le+portion03+1;
    px8=X-bar_le+portion04;
    a(py7:py8, px7:px8 )=clean_LED(py7-round(lefty_bm):py8-round(lefty_bm),  px7-round(leftx_bm):px8-round(leftx_bm)); 
%part5
    py9=Y-bar_wid;
    py10=Y+bar_wid;
    px9=X-bar_le+portion04+1;
    px10=X-bar_le+portion05;
    a(py9:py10, px9:px10 )=clean_LED(py9-round(lefty_bm):py10-round(lefty_bm),  px9-round(leftx_bm):px10-round(leftx_bm)); 
%part6
    py11=Y-bar_wid;
    py12=Y+bar_wid+1;
    px11=X-bar_le+portion05+1;
    px12=X+bar_le;
    a(py11:py12, px11:px12 )=clean_LED(py11-round(lefty_bm):py12-round(lefty_bm), px11-round(leftx_bm):px12-round(leftx_bm)); 
   
else
barX=X-round(leftx_bm);
barY=round(Y)-round(lefty_bm);
barImage=clean_LED(barY-bar_wid: barY+bar_wid, barX-bar_le: barX+bar_le); 
a(Y-bar_wid:Y+bar_wid,X-bar_le:X+bar_le)=barImage;
end


%number representing bar position
imgNum=imread(all_fileNumber(number_repeat(kk)).name); %read in image of number
imgNum2=rgb2gray(imgNum);
% figure;imshow(imgNum2)
imgNum2(imgNum2<35)=0;
imgNum2(imgNum2>=35)=255;
imgNum2(imgNum2==255)=1;
temp1=420;
a(temp1+1:temp1+43, 1232: 1280)=imgNum2;

%square_flicker
if mod(kk,3)==1 %odd number
a(500-35:500+35,1230:1280)=1; % white square
elseif mod(kk,3)==2
a(500-35:500+35,1230:1280)=0.2; %gray %test02
else
a(500-35:500+35,1230:1280)=0; % dark
end

% imwrite(a, [ori_dir, '\', num2str(kk,'%06i'),'.png']);
%imwrite conversion would change data value and class/ Indexed images are converted to RGB before writing out JPEG files

img=[];
img=a;
writeVideo(writerObj,  img);

end

%end part video
for kk=1:10
            img=0*ones(1024,1280);
            img(500-35:500+35,1230:1280)=0.2; %gray 
            writeVideo(writerObj,  img);
end
close(writerObj);

ori_dir02 =['D:\Yiko\HMM video frame\workspace_',name,'.mat'];   
save(ori_dir02)
end

 
%% HMM diagonal 19-46 (top-left/bottom right)
clear all
G_list=[ 2.5 3 4.3 4.5 5.3 6.3 6.5 7.5 9 12 20];
countt=1;
for Gvalue=G_list
clearvars -except G_list Gvalue countt
load('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\workspace_for_video')

%for randon number files
cd('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\rn_workspace') ;  
all_file = dir('*.mat') ;
file = all_file(countt).name ;
[pathstr, name, ext] = fileparts(file);
directory = [pathstr,'\'];
filename = [name,ext];
load([filename]);
name=[name];
name
Gvalue
countt=countt+1;


fps =60;  %freq of the screen flipping 
T=7*60; %second
dt=1/fps;
T=dt:dt:T;

G_HMM = Gvalue; % damping / only G will influence correlation time
D_HMM = 2700000; %dynamical range
omega =G_HMM/2.12;   % omega = G/(2w)=1.06; follow Bielak's overdamped dynamics/ 2015PNAS

%%Diagnal 19-46 > same distance as RL and UD
init=channel_map{10}+2;
final=channel_map{55}-2;

%make clean_LEDdiag : enlarge the matrix 
clean_LED=clean_LED./255; %make it to 0-1 range for double (due to imwrite format)
clean_LED(clean_LED>1)=1;

clean_LEDdiag=zeros(941,941);
clean_LEDdiag(1:200,1:200)=clean_LED(1,1);
clean_LEDdiag(201:201+length(clean_LED)-1 , 201:201+length(clean_LED)-1)=clean_LED;
for u=1:length(clean_LED)
    clean_LEDdiag(200+u,1:200)=clean_LED(u,1);%left
    clean_LEDdiag(1:200,200+u)=clean_LED(1,u); %top
    clean_LEDdiag(200+length(clean_LED)+1:400+length(clean_LED),200+u)=clean_LED(length(clean_LED),u); %bottom
    clean_LEDdiag(200+u, 200+length(clean_LED)+1:400+length(clean_LED))=clean_LED(u,length(clean_LED)); %right
end
clean_LEDdiag(201+length(clean_LED):end,1:200)=clean_LED(end,1);
clean_LEDdiag(1:200,201+length(clean_LED):end)=clean_LED(1,end);
clean_LEDdiag(201+length(clean_LED):end,201+length(clean_LED):end)=clean_LED(end,end);

%other parameters
mea_size=503;
mea_size_bm=length(clean_LEDdiag); %for clean_LEDdiag matrix
meaCenter_x=630; 
meaCenter_y=573; 
leftx_bm=meaCenter_x-(mea_size_bm-1)/2; %the first x position of the bigger mea region(luminance calibrated region) on LED screen
lefty_bm=meaCenter_y-(mea_size_bm-1)/2;

%generate location series
Yarray = zeros(1,length(T));   Xarray = zeros(1,length(T));
Vy = zeros(1,length(T));  Vx = zeros(1,length(T));
Yarray(1,1)=0;  Xarray(1,1)=0;
for t = 1:length(T)-1 
        Xarray(t+1) = Xarray(t) + Vx(t)*dt;
        Yarray(t+1) = Yarray(t) + Vy(t)*dt; 
        Vx(t+1) = (1-G_HMM*dt)*Vx(t) - omega^2*Xarray(t)*dt + sqrt(dt*D_HMM)*rntest(t); % rntest(t): gaussian noise with zero mean and unit variance
        Vy=Vx;
end
% Normalize to proper moving range
d_1=((final(1)-init(1))^2+(final(2)-init(2))^2)^(1/2);
d_2=((max(Yarray)-min(Yarray))^2+(max(Xarray)-min(Xarray))^2)^(1/2);
nr=d_1/d_2;
Xarray2=Xarray*nr;  Yarray2=Yarray*nr;
% X
mdist_x=abs(min(Xarray2)-init(1));
    if min(Xarray2)>init(1)
     mdist_x=-mdist_x;
    end
Xarray3=Xarray2+mdist_x;
newXarray=round(Xarray3);
% Y
mdist_y=abs(min(Yarray2)-init(2));
    if min(Yarray2)>init(2)
     mdist_y=-mdist_y;
    end
Yarray3=Yarray2+mdist_y;
newYarray=round(Yarray3);


%number file
number_repeat=repmat([1:12000],1,3);
cd('D:\Yiko\Timing stuff\numberFlip')
all_fileNumber= dir('*.jpeg') ;
%video frame file
name=['HMM 1946 G',num2str(G_HMM) ,' 10min Q85'];
name
% ori_dir =['D:\Yiko\HMM video frame\',name];   

%video setting
Time=T; %sec
video_fps=fps;
writerObj = VideoWriter(['D:\Yiko\HMM video frame\', name,'.avi']);  %change video name here!
writerObj.FrameRate = video_fps;
writerObj.Quality = 85;
open(writerObj);
%start part: dark adaptation
for kk=1:60
img=0*ones(1024,1280);   
writeVideo(writerObj,  img);
end


for kk =1:length(T)
a=zeros(1024,1280);%full screen pixel matrix %it's the LED screen size

%HMM 1946 bar trajectory
X=newXarray(kk);
Y=newYarray(kk);

rightopY=Y-177+8;
edge01=[];
edge01(1)=X-177-8;
edge01(2)=Y+177-8;
qoo=Y+177-8-(rightopY-1);
ptemp01=[];
for g=1:qoo
ptemp01(1)=edge01(1)+(g-1);
ptemp01(2)=edge01(2)-(g-1);
a(ptemp01(2), ptemp01(1):ptemp01(1)+32)=clean_LEDdiag(ptemp01(2)-round(lefty_bm), ptemp01(1)-round(leftx_bm):ptemp01(1)-round(leftx_bm)+32);
end

for k=1:23-1
    a(edge01(2)+k, edge01(1)+k:edge01(1)+32-k)=clean_LEDdiag(edge01(2)+k-round(lefty_bm), edge01(1)+k-round(leftx_bm):edge01(1)+32-k-round(leftx_bm)); %bottom left corner
    a(ptemp01(2)-k, ptemp01(1)+k:ptemp01(1)+32-k)=clean_LEDdiag(ptemp01(2)-k-round(lefty_bm),ptemp01(1)+k-round(leftx_bm):ptemp01(1)+32-k-round(leftx_bm)); %top right corner
end

%number representing bar position
imgNum=imread(all_fileNumber(number_repeat(kk)).name); %read in image of number
imgNum2=rgb2gray(imgNum);
% figure;imshow(imgNum2)
imgNum2(imgNum2<35)=0;
imgNum2(imgNum2>=35)=255;
imgNum2(imgNum2==255)=1;
temp1=420;
a(temp1+1:temp1+43, 1232: 1280)=imgNum2;
%square_flicker
if mod(kk,3)==1 %odd number
a(500-35:500+35,1230:1280)=1; % white square
elseif mod(kk,3)==2
a(500-35:500+35,1230:1280)=0.2; %gray %test02
else
a(500-35:500+35,1230:1280)=0; % dark
end

% imwrite(a, [ori_dir, '\', num2str(kk,'%06i'),'.png']);

img=[];
img=a;
writeVideo(writerObj,  img);
end

%end part video
for kk=1:10
            img=0*ones(1024,1280);
            img(500-35:500+35,1230:1280)=0.2; %gray 
            writeVideo(writerObj,  img);
end
close(writerObj);

ori_dir02 =['D:\Yiko\HMM video frame\workspace_',name,'.mat'];   
save(ori_dir02)
end



%% HMM 4322 diagonal
clear all
G_list=[ 2.5 3 4.3 4.5 5.3 6.3 6.5 7.5 9 12 20];
countt=1;
for Gvalue=G_list
clearvars -except G_list Gvalue countt
load('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\workspace_for_video')

%for randon number files
cd('\\192.168.0.100\Experiment\Retina\YiKo\0421 new video Br25\rn_workspace') ;  
all_file = dir('*.mat') ;
file = all_file(countt).name ;
[pathstr, name, ext] = fileparts(file);
directory = [pathstr,'\'];
filename = [name,ext];
load([filename]);
name=[name];
name
Gvalue
countt=countt+1;


fps =60;  %freq of the screen flipping 
T=7*60; %second
dt=1/fps;
T=dt:dt:T;

G_HMM =Gvalue; %Gvalue; % damping / only G will influence correlation time
D_HMM = 2700000; %dynamical range
omega =G_HMM/2.12;   % omega = G/(2w)=1.06; follow Bielak's overdamped dynamics/ 2015PNAS

%%Diagnal 4322 > same distance as RL and UD
init(1)=channel_map{50}(1)-2;
init(2)=channel_map{50}(2)+2;
final(1)=channel_map{15}(1)+2;
final(2)=channel_map{15}(2)-2;

%make clean_LEDdiag : enlarge the matrix 
clean_LED=clean_LED./255; %make it to 0-1 range for double (due to imwrite format)
clean_LED(clean_LED>1)=1;

clean_LEDdiag=zeros(941,941);
clean_LEDdiag(1:200,1:200)=clean_LED(1,1);
clean_LEDdiag(201:201+length(clean_LED)-1 , 201:201+length(clean_LED)-1)=clean_LED;
for u=1:length(clean_LED)
    clean_LEDdiag(200+u,1:200)=clean_LED(u,1);%left
    clean_LEDdiag(1:200,200+u)=clean_LED(1,u); %top
    clean_LEDdiag(200+length(clean_LED)+1:400+length(clean_LED),200+u)=clean_LED(length(clean_LED),u); %bottom
    clean_LEDdiag(200+u, 200+length(clean_LED)+1:400+length(clean_LED))=clean_LED(u,length(clean_LED)); %right
end
clean_LEDdiag(201+length(clean_LED):end,1:200)=clean_LED(end,1);
clean_LEDdiag(1:200,201+length(clean_LED):end)=clean_LED(1,end);
clean_LEDdiag(201+length(clean_LED):end,201+length(clean_LED):end)=clean_LED(end,end);

%other parameters
mea_size=503;
mea_size_bm=length(clean_LEDdiag); %for clean_LEDdiag matrix
meaCenter_x=630; 
meaCenter_y=573; 
leftx_bm=meaCenter_x-(mea_size_bm-1)/2; %the first x position of the bigger mea region(luminance calibrated region) on LED screen
lefty_bm=meaCenter_y-(mea_size_bm-1)/2;

%generate location series
Yarray = zeros(1,length(T));   Xarray = zeros(1,length(T));
Vy = zeros(1,length(T));  Vx = zeros(1,length(T));
Yarray(1,1)=0;  Xarray(1,1)=0;
for t = 1:length(T)-1 
        Xarray(t+1) = Xarray(t) + Vx(t)*dt;
        Yarray(t+1) = Yarray(t) + Vy(t)*dt; 
        Vx(t+1) = (1-G_HMM*dt)*Vx(t) - omega^2*Xarray(t)*dt + sqrt(dt*D_HMM)*rntest(t); % rntest(t): gaussian noise with zero mean and unit variance
        Vy=-Vx;
end
% Normalize to proper moving range
d_1=((final(1)-init(1))^2+(final(2)-init(2))^2)^(1/2);
d_2=((max(Yarray)-min(Yarray))^2+(max(Xarray)-min(Xarray))^2)^(1/2);
nr=d_1/d_2;
Xarray2=Xarray*nr;  Yarray2=Yarray*nr;
% X
mdist_x=abs(min(Xarray2)-final(1));
    if max(Xarray2)>init(1)
     mdist_x=-mdist_x;
    end
Xarray3=Xarray2+mdist_x;
newXarray=round(Xarray3);
% Y
mdist_y=abs(min(Yarray2)-init(2));
    if min(Yarray2)>init(2)
     mdist_y=-mdist_y;
    end
Yarray3=Yarray2+mdist_y;
newYarray=round(Yarray3);

%number file
number_repeat=repmat([1:12000],1,3);
cd('D:\Yiko\Timing stuff\numberFlip')
all_fileNumber= dir('*.jpeg') ;
%video frame file
name=['HMM 4322 G',num2str(G_HMM) ,' 10min Q85'];
name
% ori_dir =['D:\Yiko\HMM video frame\',name];   

%video setting
Time=T; %sec
video_fps=fps;
writerObj = VideoWriter(['D:\Yiko\HMM video frame\', name,'.avi']);  %change video name here!
writerObj.FrameRate = video_fps;
writerObj.Quality = 85;
open(writerObj);
%start part: dark adaptation
for kk=1:60
img=0*ones(1024,1280);   
writeVideo(writerObj,  img);
end



for kk =1: length(T)
a=zeros(1024,1280);%full screen pixel matrix %it's the LED screen size

%HMM 4322 bar trajectory
X=newXarray(kk);
Y=newYarray(kk);

rightbY=Y+177-8;
edge01=[];
edge01(1)=X-177-8;
edge01(2)=Y-177+8;
qoo=rightbY-(edge01(2)-1);
ptemp01=[];
for g=1:qoo
ptemp01(1)=edge01(1)+(g-1);
ptemp01(2)=edge01(2)+(g-1);
a(ptemp01(2), ptemp01(1):ptemp01(1)+32)=clean_LEDdiag(ptemp01(2)-round(lefty_bm), ptemp01(1)-round(leftx_bm):ptemp01(1)-round(leftx_bm)+32);
end

for k=1:23-1
    a(edge01(2)-k, edge01(1)+k:edge01(1)+32-k)=clean_LEDdiag(edge01(2)-k-round(lefty_bm), edge01(1)+k-round(leftx_bm):edge01(1)+32-k-round(leftx_bm)); %bottom left corner
    a(ptemp01(2)+k, ptemp01(1)+k:ptemp01(1)+32-k)=clean_LEDdiag(ptemp01(2)+k-round(lefty_bm),ptemp01(1)+k-round(leftx_bm):ptemp01(1)+32-k-round(leftx_bm)); %top right corner
end

%number representing bar position
imgNum=imread(all_fileNumber(number_repeat(kk)).name); %read in image of number
imgNum2=rgb2gray(imgNum);
% figure;imshow(imgNum2)
imgNum2(imgNum2<35)=0;
imgNum2(imgNum2>=35)=255;
imgNum2(imgNum2==255)=1;
temp1=420;
a(temp1+1:temp1+43, 1232: 1280)=imgNum2;

%square_flicker
if mod(kk,3)==1 %odd number
a(500-35:500+35,1230:1280)=1; % white square
elseif mod(kk,3)==2
a(500-35:500+35,1230:1280)=0.2; %gray %test02
else
a(500-35:500+35,1230:1280)=0; % dark
end

% imwrite(a, [ori_dir, '\', num2str(kk,'%06i'),'.png']);
img=[];
img=a;
writeVideo(writerObj,  img);

end

%end part video
for kk=1:10
            img=0*ones(1024,1280);
            img(500-35:500+35,1230:1280)=0.2; %gray 
            writeVideo(writerObj,  img);
end
close(writerObj);

ori_dir02 =['D:\Yiko\HMM video frame\workspace_',name,'.mat'];   
save(ori_dir02)
end









% %% Generate video
% clear all
% Time=5*60; %sec
% video_fps=60;
% 
% writerObj = VideoWriter('D:\Yiko\HMM video frame\HMM RL G5 5min02.avi');  %change video name here!
% writerObj.FrameRate = video_fps;
% open(writerObj);
% 
% %start part: dark adaptation
% for kk=1:60
% img=0*ones(1024,1280);   
% writeVideo(writerObj,  img);
% end
% 
% %stimulation frame
% cd('D:\Yiko\HMM video frame\HMM RL G5 5min02') %load in the image files
% all_file = dir('*.jpeg') ;
% for i=1:Time*video_fps
%             img=imread(all_file(i).name);
%             writeVideo(writerObj,  img);
% end
% 
% %end part
% for kk=1:5
%             img=0.2*ones(1024,1280); 
%             writeVideo(writerObj,  img);
% end
% close(writerObj);


%% Ref code
% X=newXarray(j);
% a(Y-bar_le:Y+bar_le,X-bar_wid:X+bar_wid)=255;
% 
% barX=newXarray(1,runi)-round(leftx_bm);
% barY=round(Y)-round(lefty_bm);
% barImage=clean_LED(barY-bar_le: barY+bar_le, barX-bar_wid: barX+bar_wid); 

