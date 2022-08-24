%ImExtrapolation_function
%%Quadratic extrapolation for finite data to estimate the real mutual infromation value; input includes a 2*N data set and the states for partition
%%free parameters: Ns fractions selected; rep repitition for smoothing data; rescalong the fitting function or deviding the time window beforehand

function ImEx = ImExtrapolation_function(dat,states)

rep = 30;% repeat for every point
ImEx = zeros(1,rep);
for rr = 1:rep
    Ns = 5; frac = 1:Ns;
    fraction = 1./frac;
    A = [ones(1,length(frac));  1./fliplr(frac);  1./fliplr(frac).^2];  %*.1*.01
    naive = zeros(1,Ns);
    
    for f = 1:length(fraction)
        sampN = floor(size(dat,1)*fraction(f));
        p = randperm(size(dat,1));
        subsample = dat(p(1:sampN),:);
        
        [N,C] = hist3(subsample,[states,states]);    
        px = sum(N,1)/sum(sum(N));
        py = sum(N,2)/sum(sum(N));
        pxy = N/sum(sum(N));
        temp = [];
        for j = 1:length(px)
            for k = 1:length(py)
                temp(k,j) = pxy(k,j)*log( pxy(k,j)/ (py(k)*px(j)) )/log(2);  %/(bin/1000);
            end
        end
        naive(f) = nansum(temp(:));
    end
    [pp,S] = polyfit(frac,naive,2);
    fitting = pp(end);  %A'\naive';  %
    ImEx(rr) = fitting(1);%max(0,fitting(1));  %
end
ImEx = mean(ImEx);

end


