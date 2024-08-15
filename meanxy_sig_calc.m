function [meanxy_sig] = meanxy_sig_calc(dmean)
if dmean == 0
    meanxy_sig=0;
else
    meanxy_sig=sqrt(sum(dmean.^2)./(length(dmean)));
end