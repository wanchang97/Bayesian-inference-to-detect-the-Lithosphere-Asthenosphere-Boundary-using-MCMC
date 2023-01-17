function M = ranodmWalkSampler(prior,Likelihood,m0,nRW,sigmaRW_ratio)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generate samples following of a Markov chain whose stationary distribution
is the posterior distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    prior: 
    likelihood:
    m0              : a point chosen deterministically
    nRW: 
    sigmaRW_ratio   : sqrt of diagonal of the covariance matrix|dim :np x 1
Outputs: 
    M               :  Samples of the markov chain | dim: np x nsamples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%}

% Idea 1: specify the number of samples nRW 
np = size(m0,1);
M = zeros(np,nRW);
M(:,1) = m0;
muRW = zeros(np,1);
sigmaRW = m0 *sigmaRW_ratio;
covRW = diag(sigmaRW).^2;
for i = 2 : nRW
    % generate the candidate sample based on the given sample
    z = mvnrnd(muRW,covRW,1)';
    x = M(:,i-1) + z;
    % convergence ratio
    r = prior(x)*Likelihood(x)/(prior(M(:,i-1))*Likelihood(M(:,i-1)));
    alpha = min(1,r);
    % generate a sample u from the uniform distribution
    u = unifrnd(0,1,1);
    if u <= alpha % accept the sample
        M(:,i) = x;
    else % reject the sample
        M(:,i) = M(:,i-1);
    end
end
% the burn-in period is set to be 0.1 *nRW
nburn = nRW * 0.1;
M = M(:,nburn:nRW);

% Idea 2: do not force how many samples to generate, but use a tolerance to
% control the convergence
