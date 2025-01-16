function [ms, vs, ss, mds,xscale] = movstats(x,w,ws),

x  = x(:);
ii = 0;
for i=1:ws:length(x);
    ii = ii + 1;
    if i<w/2
        x1=x(1:i+w/2);
    elseif abs(length(x)-i)<w/2
        x1=x(i-w/2+1:end);
    else
        x1=x(i-w/2+1:i+w/2);
    end
    vs(ii,1) = var(x1);
    ms(ii,1) =  mean(x1);
    mds(ii,1) =  median(x1);
    ss(ii,1) = std(x1);
    xscale(ii,1) = i;
end

