% % Confidence interval for unknown mean and known standard deviation
% samples = [102.5,101.7,103.1,100.9,100.5,102.2];
% samplemean = mean(samples)
% samplestd = std(samples)
% samplestd2 = sqrt(sum((samples-samplemean).^2)/(length(samples)-1))
% samplestd3 = sqrt(sum((samples-samplemean).^2)/(length(samples)))
% % Confidence interval for unknown mean and unknown standard deviation
function delta = CIinterval(alpha,samples)
% This function is to calculate the confidence interval of a group of
% samples with unknown mean and unknown standard deviations. Note that due
% to the central limit theorem, we assume that the samples follow the
% normal distribution (z- distribution)
% Input
%   alpha : 0.05 -- 95% CI
%   samples : dim nsample x np
% Output
%   delta: dim np x 1
%
sample_std = std(samples,1);
nsamples = size(samples,1);
zstar = norminv(1-alpha/2);
zstar2 = zstar *2.3;
zstar3 = tinv(1-alpha/2,nsamples);
%delta = sample_std/sqrt(nsamples)*zstar;
%delta = sample_std*zstar;
%delta = sample_std/sqrt(nsamples)*zstar2;
delta = sample_std/sqrt(nsamples)*zstar3;
delta = delta';

end