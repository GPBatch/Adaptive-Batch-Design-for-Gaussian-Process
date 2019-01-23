function [lee, er, ee, bias] = gp_perf(Ef, Varf, xt, ft, pcr, lambda)

%%%%%%%%%%% This function is used to calculate the performance metric for
%%%%%%%%%%% level set estimation. (local/integral empirical error and bias)
d = size(xt,2);
m = size(xt,1);
nt1 = lambda*m;

lee = metricmee(Ef, Varf);

if ( d == 1 )
    ee = sum(lee)/m;
    er = sum(Ef.*ft<0)/m;
    bias = (sum((Ef>0)&(ft<0))-sum((Ef<0)&(ft>0)))/m;
else
    er = pcr*sum(Ef(1:nt1).*ft(1:nt1)<0)/nt1 + (1-pcr)*sum(Ef((nt1+1):end).*ft((nt1+1):end)<0)/(m - nt1);
    ee = pcr*sum(lee(1:nt1))/nt1+(1-pcr)*sum(lee((nt1+1):end))/(m - nt1);
    bias = pcr*(sum((Ef(1:nt1)>0)&(ft(1:nt1)<0))-sum((Ef(1:nt1)<0)&(ft(1:nt1)>0)))/nt1 + (1-pcr)*(sum((Ef((nt1+1):end)>0)&(ft((nt1+1):end)<0))-sum((Ef((nt1+1):end)<0)&(ft((nt1+1):end)>0)))/(m - nt1);
end

end