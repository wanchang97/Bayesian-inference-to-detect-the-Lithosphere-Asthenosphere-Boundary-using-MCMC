function [M,acc] = logAM(logpiPdf,loglikelihood,m0,nRW,sigmam,fwd)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generate samples following of a Markov chain whose stationary distribution
is the posterior distribution
Adaptive Metropolis algorithm
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
    acc             :The number of samples being accepted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
k = fwd.k;
nsd = fwd.nsd;
switch fwd.mtype
    case 0
        d0 = m0;
        sigmad = sigmam.d;
        D = zeros(nRW,k,1);
        D(1,:) = d0;
        muDRW = zeros(k,1);
        covDRW = diag(sigmad.^2*ones(1,k));
        acc = 0;
        tune_d = 2.38^2/k;
        for i = 2:nRW
            if rem(i,10) == 1
                fprintf('Regular grids RW sample: %d \n ', i )
            end
            % update the proposal covariance matrix
            if i > 2 % update the covariance matrix by
                covD_emp = empiricalCovariance(D(1:i-1,:));
                covDRW = tune_d*covD_emp;
                % to avoid sigularity
                covDRW  = covDRW + diag(10^(-4)*ones(1,k));
            end
            %propose a new candidate from the proposal density gaussian
            z_depth = mvnrnd(muDRW,covDRW)';
            d_candidate = D(i-1,:)' + z_depth;
            r = loglikelihood(d_candidate)+logpiPdf(d_candidate)-(loglikelihood(D(i-1,:)')+logpiPdf(D(i-1,:)'));
            alpha = min(0,r);
            % generate a sample u from the uniform distribution
            u = unifrnd(0,1,1);
            u = log(u);
            if u <= alpha % accept the sample
                D(i,:) = d_candidate;
                acc = acc + 1;
            else % reject the sample
                D(i,:) = D(i-1,:);
            end
            M = D;
        end
    case 1
        b0 = m0(1:k-1);
        d0 = m0(k:end);
        sigmab = sigmam.b;
        sigmad = sigmam.d;
        B = zeros(nRW,k-1,nsd-1);
        D = zeros(nRW,k,1);
        B(1,:) = b0;
        D(1,:) = d0;

        muBRW = zeros((k-1)*(nsd-1),1);
        covBRW = diag(sigmab.^2*ones(1,k-1));
        muDRW = zeros(k,1);
        covDRW = diag(sigmad.^2*ones(1,k));

        acc = 0;
        tune_b = 2.38^2/(k-1);
        tune_d = 2.38^2/k;
        for i = 2:nRW
            if rem(i,10) == 1
                fprintf('Irregular grids RW sample: %d \n ', i )
            end
            % update the proposal covariance matrix
            if i > 2 % update the covariance matrix by
                covB_emp = empiricalCovariance(B(1:i-1,:));
                covBRW = tune_b*covB_emp;
                covD_emp = empiricalCovariance(D(1:i-1,:));
                covDRW = tune_d*covD_emp;
                % to aboid sigularity
                covBRW  = covBRW + diag(10^(-4)*ones(1,k-1));
                covDRW  = covDRW + diag(10^(-4)*ones(1,k));
            end
            %propose a new candidate from the proposal density gaussian
            z_b = mvnrnd(muBRW,covBRW)';
            z_depth = mvnrnd(muDRW,covDRW)';
            b_candidate = B(i-1,:)' + z_b;
            d_candidate = D(i-1,:)' + z_depth;
            r = loglikelihood([b_candidate;d_candidate])+logpiPdf([b_candidate;d_candidate])-(loglikelihood([B(i-1,:)';D(i-1,:)'])+logpiPdf([B(i-1,:)';D(i-1,:)']));
            alpha = min(0,r);
            % generate a sample u from the uniform distribution
            u = unifrnd(0,1,1);
            u = log(u);
            if u <= alpha % accept the sample
                B(i,:) = b_candidate;
                D(i,:) = d_candidate;
                acc = acc + 1;
            else % reject the sample
                B(i,:) = B(i-1,:);
                D(i,:) = D(i-1,:);
            end
            M =[B,D];
        end

end





