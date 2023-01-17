function [samplesM,Z] = modifiedTMCMC(priorInfo,fwd,likelihood,n,s)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modified TMCMC, using the resampling step to sample the intermediate samples
in linear scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    priorInfo       : has the distribution of inputs
    fwd             :   
    likelihood:     :
    n               :The number of samples generated in each step, totoally nl samples
    s               :The total number of intermediate levels
    sigmad          :sqrt of diagonal of the covariance matrix|dim :k x 1
Outputs: 
    samplesM        :Samples of all s stages |dim: s x n x nm
    Z               :Evidence estimation |dim: 1x1
    acc             :The number of samples being accepted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%% Some contants and initialization
beta        = 0.2;          % Prescribed scaling factor (0.2 recommended choice)
thres_cov   = 1;            % threshold for the c.o.v (100% recommended choice)
% memory allocation
samplesM    = cell(s,1);    % spaces for the samples
Z_t         = ones(s,1);    % spaces for intermediate evidence
p           = zeros(s,1);   % store tempering parameters

% 1. Obtain n samples from the prior PDF and evaluate likelihood
m0 = random(priorInfo.mdist,1,n);
samplesM{1} = m0;
w = zeros(s,n); % Initialize the weights
L_t = zeros(n,1);
for i= 1:n
    L_t(i) = likelihood(m0(1,i));
end
%% TMCMC
for t  =2:s
    fprintf('\n regular grids modified TMCMC intermediate level t = %g, with p_{t} = %g\n', t-1, p(t-1));
    % 2. Compute tempering parameter p_{t}
    % dp = p_{t}-p_{t-1}
    fun = @(dp) std(L_t.^dp)/mean(L_t.^dp) -thres_cov;
    [dp,~,flag] = fzero(fun,0);
    % if fzero doesn't work, try with fsolve
    if flag > 0 % work
        dp = abs(dp); % dp is >=0
    elseif license('test','optimization_toolbox')
        option = optimset('Display','off');
        [dp,~,flag] = fsolve(fun,0,option);
        dp = abs(dp); % dp is <=0
        if flag < 0
            error('fzero and fsolve do not converge');
        end
    else error('no optimization_toolbox available');
    end
    if ~isnan(dq) % when dq is a finite value
        p(t) = min(1,p(t-1)+dp);
    else % when dq is NaN infinite value
        p(t) = 1;
        fprintf('Variable p was set to %f, since it is not possible to find a suitable value \n',p(t));
    end
    % 3. Compute the weights w_t and evidence Z_t
    w_t = L_t^(p(t)-p(t-1));
    Z_t(t) = mean(w_t);
    w_t_norm = w_t/sum(w_t);
    weightscumsum = cumsum(w_t_norm);
    % 4 Calculating the estimated covariance matrix of f_t(m) Gaussian
    % proposal PDF based on samples samplesM{t-1}

    temp = w_t * samplesM{t-1};
    temp2 = (samplesM{t-1}-temp/sum(w_t)) * (samplesM{t-1}-temp/sum(w_t))';
    for i = 1:n % n_{t-1}
        weightedsum = weightedsum +  w_t(i) * temp2(i);% ?
    end
    sigma_estimated = beta^2 * weightedsum;
    pdf = % Gaussian proposal pdf
    proprnd =  % Covariance matrix equal to sigma_estimated
    % 5 Metropolis resampling step to obtain n samples following f_t(m)
    
    X_update = samplesM{t-1};
    index = zeros(n,1);
    for i = 1:n
        index(i) = find(rand <= weightscumsum,1);
        x_chosen = X_update(index(i),:);
        x_MHsample = mhsample(x_chosen,1,'pdf',pdf,'proprnd',proprnd,'symmetric',1);
        X_update(i,1) = x_MHsample;
        X_chosen(i,:) = x_chosen;
        samplesM{t}(i) = x_MHsample;
        L_t(i) = likelihood(samplesM{t}(1,i));
    end
end
Z = prod(Z_t);
