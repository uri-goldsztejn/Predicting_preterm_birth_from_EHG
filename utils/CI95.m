function [CI] = CI95(data)
% CI95 Computes the 95% confidence interval of the mean of the data. 

SEM = std(data)/sqrt(length(data));               % Standard Error
ts = tinv([0.025  0.975],length(data)-1);      % T-Score
CI = mean(data) + ts*SEM;
end