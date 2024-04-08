function sample=emergencesample(T)
% we take advantage of the fact that the cumulative distribution of the emergence time is exp(t)/(exp(t)+exp(T)), 
% hence the t that corresponds to a uniformly distributed p is ln(exp(T) p/(1-p)) = ln(exp(T))+ln(p/(1-p)) = T+ln(p/(1-p))

p=rand(size(T));
sample=T+log(p./(1-p));

