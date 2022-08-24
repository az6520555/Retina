function[t,information,corr] = MIfun(Neurons,isi2,bin,backward,forward);
%% variables
% bin=1; % ms  
BinningInterval = bin*10^-3;
% backward=ceil(1000/bin); forward=ceil(1000/bin); 

dat=[];informationp=[];temp=backward+2;
    for i=1:backward+1 %past(t<0)
        x = Neurons((i-1)+forward+1:length(Neurons)-backward+(i-1))';
        y = isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        norm=1;

        [N,C]=hist3(dat{i}); 
        px=sum(N,1)/sum(sum(N)); 
        py=sum(N,2)/sum(sum(N)); 
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
              temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2)/norm;
            end
        end
        temp=temp-1;
        informationp(temp)=nansum(temp2(:));
        c=corrcoef(x,y);
        corrp(temp)=c(2,1);
    end  

    dat=[];informationf=[];temp=0;sdat=[];
    for i=1:forward
        x =Neurons(forward+1-i:length(Neurons)-backward-i)';
        y = isi2(forward+1:length(isi2)-backward)';
        dat{i}=[x,y];
        norm=1;

        [N,C]=hist3(dat{i}); 
        px=sum(N,1)/sum(sum(N)); 
        py=sum(N,2)/sum(sum(N)); 
        pxy=N/sum(sum(N));
        temp2=[];
        for j=1:length(px)
            for k=1:length(py)
                temp2(k,j)=pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2)/norm;
            end
        end
        temp=temp+1;
        informationf(temp)=nansum(temp2(:));
        c=corrcoef(x,y);
        corrf(temp)=c(2,1);
    end
%     information=[informationp informationf]/BinningInterval;
    information=[informationp informationf];
    corr = [corrp corrf];
    t=[-backward*bin:bin:forward*bin];  
%     figure;hold on;plot(t,information,'LineWidth',2);
%     xlabel('\deltat');ylabel('MI (bits/s)');