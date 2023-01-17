function [M,lnL,acc] = logRW(logpiPdf,loglikelihood,m0,nRW,sigmam,fwd)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generate samples following of a Markov chain whose stationary distribution
is the posterior distribution
Metropolis algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    logpiPdf        : 
    loglikelihood:  :
    c0              :The initial sites location vector |dim k x nsd
    d0              :The initial LAB depth vector |dim k-1 x 1
    nRW             :The number of maximal interations
    sigmab          :sqrt of diagonal of the covariance matrix|dim :(k-1) x 1
    sigmad          :sqrt of diagonal of the covariance matrix|dim :k x 1
Outputs: 
    B               :Samples of sites in Markov chain |dim: nsamples x (k x nsd)
    D               :Samples of depth in Markov chain |dim: nsamples x (k x 1)
    lnL             :evaluated Loglikelihood for each accepted samples
    acc             :The number of samples being accepted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
k = fwd.k;
nsd = fwd.nsd;
lnL = zeros(nRW,1);
D = zeros(nRW,k,1);
switch fwd.mtype
    case 0
        d0 = m0;
        sigmad = sigmam.d;  
        D(1,:) = d0;
        lnL(1) = loglikelihood(D(1,:)');
        muDRW = zeros(k,1);
        covDRW = diag(sigmad.^2*ones(1,k));
        acc = 0;
        for i = 2:nRW
            if rem(i,10) == 1
                fprintf('Regular grids RW sample: %d \n ', i )
            end
            % propose a new candidate coordinate from the proposal density
            z_depth = mvnrnd(muDRW,covDRW)';
            d_candidate = D(i-1,:)' + z_depth;
            loglikelihood_candidate = loglikelihood(d_candidate);
            r = loglikelihood_candidate+logpiPdf(d_candidate)-(lnL(i-1)+logpiPdf(D(i-1,:)'));
            alpha = min(0,r);
            % generate a sample u from the uniform distribution
            u = unifrnd(0,1,1);
            u = log(u);
            if u <= alpha % accept the sample
                D(i,:) = d_candidate;
                acc = acc + 1;
                lnL(i) = loglikelihood_candidate;
            else % reject the sample
                D(i,:) = D(i-1,:);
                lnL(i) = lnL(i-1);
            end
            M = D;
        end
    case 1
        b0 = m0(1:k-1);
        d0 = m0(k:end);
        sigmab = sigmam.b;
        sigmad = sigmam.d;
        B = zeros(nRW,k-1,nsd-1);
        B(1,:) = b0;
        D(1,:) = d0;
        lnL(1) = loglikelihood([B(1,:)';D(1,:)']);
        muBRW = zeros((k-1)*(nsd-1),1);
        covBRW = diag(sigmab.^2*ones(1,k-1));
        muDRW = zeros(k,1);
        covDRW = diag(sigmad.^2*ones(1,k));
        acc = 0;
        for i = 2:nRW
            if rem(i,10) == 1
                fprintf('Irregular grids RW sample: %d \n ', i )
            end
            % propose a new candidate coordinate from the proposal density
            z_coordinate = mvnrnd(muBRW,covBRW)';
            % How to make the new candidate inside the domain?
            b_candidate = B(i-1,:)' + z_coordinate;
            % and update its depth
            z_depth = mvnrnd(muDRW,covDRW)';
            %sz_depth = mvnrnd(muDRW,covDRW,1,1);
            d_candidate = D(i-1,:)' + z_depth;
            % solve the forward problem and the convergence ratio
            %r = prior(c_candidate,d_candidate)*likelihood(c_candidate,d_candidate)/(prior(C(i-1,:,:),D(i-1,:,:))*likelihood(C(i-1,:,:),D(i-1,:,:)));
            % for the uniform prior distribution
            loglikelihood_candidate = loglikelihood([b_candidate;d_candidate]);
            r = loglikelihood_candidate+logpiPdf([b_candidate;d_candidate])-(lnL(i-1)+logpiPdf([B(i-1,:)';D(i-1,:)']));
            alpha = min(0,r);
            % generate a sample u from the uniform distribution
            u = unifrnd(0,1,1);
            u = log(u); % must be negative
            if u <= alpha % accept the sample
                B(i,:) = b_candidate;
                D(i,:) = d_candidate;
                acc = acc + 1;
                lnL(i) = loglikelihood_candidate;
            else % reject the sample
                B(i,:) = B(i-1,:);
                D(i,:) = D(i-1,:);
                lnL(i) = lnL(i-1);
            end
            M = [B,D];
        end
end
% the burn-in period is set to be 0.1 *nRW
% nburn = nRW * 0.1;
% B = B(nburn:nRW,:,:);
% D = D(nburn:nRW,:,:);
end
