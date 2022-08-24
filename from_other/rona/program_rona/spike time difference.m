load('E:\google_rona\20161122\HMM3_diffr\HMM3_r=0.10_20161122.mat');

N = Spikes{6};
P = Spikes{38};

dif=[];temp = 0;
for i = 1:length(N)
    for j = 1:length(P)
        if abs(N(i)-P(j))<0.05
            temp = temp+1;
            dif(temp) = P(j)-N(i);
        end
    end
end

figure;hist(dif,30)
mean(dif)