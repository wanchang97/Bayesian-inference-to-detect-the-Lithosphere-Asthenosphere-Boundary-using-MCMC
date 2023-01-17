function [M,E,Cov] = momentum_by_RW(prior,likelihood,m0,nRW,sigmaRW_ratio,histplot,figurepath)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Evalute the quantity of interest using RW sampler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    prior: 
    likelihood:
    m0              : a point chosen deterministically
    nRW: 
    sigmaRW         : sqrt of diagonal of the covariance matrix|dim :np x 1
    histplot        : 1 plot, 0 no plot
    figurepath      : the location to save the file
Outputs: 
    E               : Expectation | dim: np x 1
    Cov             : Covariance matrix | dim: np x np
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%}
% Stage 1: Generate the samples following Markov Chain

M = ranodmWalkSampler(prior,likelihood,m0,nRW,sigmaRW_ratio); % np x (nRW-nburn)

% Stage 2: Use the stationary samples to evaluate the interested value
% expectation or variance
E = mean(M')'; % equivalent to mean(M,2)
%s = std(M')';
Cov = (M-E)*(M-E)'/nRW;
% histogram plot
if histplot == 1
    histplot_results_RW(M,figurepath,nRW,sigmaRW_ratio);
else
    return
end
