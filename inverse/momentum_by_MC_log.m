function [E,Cov,E_uncertainty] = momentum_by_MC_log(k,nMC,priorInfo,my_loglikelihood,fwd)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Estimating the expectation and variance using the crude MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    k:              number of cells
    nMC:            number of samples to be evaluated in crude MC
    priorInfo:      containing information of priorPdf
    my_likelihood: 
Outputs: 
    E:               Expectation| dim: np x nMC
    V:               Variance | dim: np x np x nMC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
% generate C samples acoording to the prior distribution
switch fwd.mtype
    case "0"% 0: Fixed discretization,fixed dimension m = d : np = k
    case "1"% 1: Updated discretization, fixed dimension m= (C,d): np = k +(nsd-1)*(k-1);

        np = k + (fwd.nsd-1)*(k-1);
        sampled_B = priorInfo.bdist.random(nMC,k-1);
        sampled_D = priorInfo.ddist.random(nMC,k);

        % evaluate the likelihood, mLikelihood at the generated samples
        sum_leval = 0;
        sum_mleval = 0;
        Leval  = zeros(nMC,1);
        for i = 1:nMC
            if rem(i,10) ==1
                fprintf("MC sample: %d \n",i);
            end
            logleval = my_loglikelihood(sampled_B(i,:)',sampled_D(i,:)');
            Leval(i,:) = logleval;
            sum_leval = sum_leval + leval;
            sum_mleval = sum_mleval + [sampled_B(i,:)';sampled_D(i,:)']*leval;
        end
        E  = sum_mleval /sum_leval;
        mmLeval = zeros(nMC,k+(k-1)*(fwd.nsd-1),k+(k-1)*(fwd.nsd-1));
        for I = 1:nMC
            mmLeval(I,:,:) = ([sampled_B(I,:)';sampled_D(I,:)']-E)*([sampled_B(I,:)';sampled_D(I,:)']-E)'*Leval(I,1);
        end
        Cov  = squeeze(sum(mmLeval,1)/sum(Leval,1));
        E_uncertainty = sqrt(diag(Cov)/nMC);
    case "2"% 2: Updated discretization,fixed dimenison m = (C,d): np = k + (nsd-1)*k;
        np = k + (fwd.nsd-1)*k;
        sampled_C = priorInfo.cdist.random(nMC,k);
        sampled_D = priorInfo.ddist.random(nMC,k);

        % evaluate the likelihood, mLikelihood at the generated samples
        sum_leval = 0;
        sum_mleval = 0;
        Leval  = zeros(nMC,1);
        sampled_M = zeros(nMC,np);
        for i = 1:nMC
            if rem(i,10) ==1
                fprintf("MC sample: %d \n",i);
            end
            leval = my_likelihood(sampled_C(i,:),sampled_D(i,:));
            Leval(i,:) = leval;
            sum_leval = sum_leval + leval;
            sampled_M(i,:) = [sampled_C(i,:)';sampled_D(i,:)'];
            sum_mleval = sum_mleval + [sampled_C(i,:)';sampled_D(i,:)']*leval;
        end
        E = sum_mleval /sum_leval;
        E2 = sampled_M*Leval/sum(Leval,1);
        % How to calculate the uncertainty?
        mmLeval = zeros(nMC,k*fwd.nsd,k*fwd.nsd);
        for I = 1:nMC
            mmLeval(I,:,:) = ([sampled_C(I,:)';sampled_D(I,:)']-E)*([sampled_C(I,:)';sampled_D(I,:)']-E)'*Leval(I,1);
        end
        Cov  = squeeze(sum(mmLeval,1)/sum(Leval,1));
    case "3"% 3: Updated discretization, fixed dimension m = (C,d):fix only the first voronoi cell: np = k +(nsd-1)*(k-1)
        np = k + (fwd.nsd-1)*(k-1);
        sampled_C = priorInfo.cdist.random(nMC,k-1);
        sampled_D = priorInfo.ddist.random(nMC,k);

        % evaluate the likelihood, mLikelihood at the generated samples
        sum_leval = 0;
        sum_mleval = 0;
        Leval  = zeros(nMC,1);
        for i = 1:nMC
            if rem(i,10) ==1
                fprintf("MC sample: %d \n",i);
            end
            leval = my_likelihood([fwd.c1;sampled_C(i,:)'],sampled_D(i,:)');
            Leval(i,:) = leval;
            sum_leval = sum_leval + leval;
            sum_mleval = sum_mleval + [sampled_C(i,:)';sampled_D(i,:)']*leval;
        end
        E  = sum_mleval /sum_leval;
        mmLeval = zeros(nMC,k+(k-1)*(fwd.nsd-1),k+(k-1)*(fwd.nsd-1));
        for I = 1:nMC
            mmLeval(I,:,:) = ([sampled_C(I,:)';sampled_D(I,:)']-E)*([sampled_C(I,:)';sampled_D(I,:)']-E)'*Leval(I,1);
        end
        Cov  = squeeze(sum(mmLeval,1)/sum(Leval,1));
    case "4"% 4: Updated discretization, transdimension m = (k,C,d): np = 1+ k + (nsd-1)*k;
    case "5"% 5: Updated discretization, transdimenison,fix the first voronoi cell m = (k,C,d): np = 1+ k + (nsd-1)*(k-1);
end