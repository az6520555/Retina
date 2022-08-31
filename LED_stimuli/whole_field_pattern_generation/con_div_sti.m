% diverge and converge stimulus
clear
sti_type={'on','div'}; % on or off, div or conv
image_time=0.1;
Ttot=300;
frame_rate=1/60;
frame_per_image=image_time/frame_rate;
nrow=Ttot/image_time;
ncol=20;
OLED_res=189;
pixel_length=int32(OLED_res/ncol);
x=zeros(nrow,ncol);
x(1,:)=randi([0,1],[1,ncol]);
x(2:end,1)=randi([0,1],[nrow-1,1]);

for i=1:ncol-1
    for j=1:nrow-1
        if strcmp(sti_type{1},'on')
            if x(j,i)+x(j,i+1)==2 || x(j,i)+x(j,i+1)==0
                x(j+1,i+1)=1;
            else 
                x(j+1,i+1)=0;
            end
        elseif strcmp(sti_type{1},'off')
            if x(j,i)+x(j,i+1)==2 || x(j,i)+x(j,i+1)==0
                x(j+1,i+1)=0;
            else 
                x(j+1,i+1)=1;
            end
        end
    end
end
%% zoom x dimension
x_zoom=zeros(nrow,ncol*pixel_length);
for i=1:ncol
    for num_pixel=1:pixel_length
        i_col=(i-1)*pixel_length+num_pixel;
        x_zoom(:,i_col)=x(:,i);
    end
end
%% zoom t dimension
x_zoom2=zeros(nrow*frame_per_image,ncol*pixel_length);
for i=1:nrow
    for j=1:frame_per_image
        i_row=(i-1)*frame_per_image+j;
        x_zoom2(i_row,:)=x_zoom(i,:);
    end
end
if strcmp(sti_type{2},'div')
    x_zoom2=flipud(x_zoom2); 
end
x_zoom2(:,OLED_res+1:end)=[];
imshow(x_zoom2')

% save(['F:\§Úªº¶³ºÝµwºÐ\Retina exp\2021 research assistant\spatial_sitmuli\new\',sti_type{2},'_',sti_type{1},'_',num2str(image_time),'.mat'],'x_zoom2')