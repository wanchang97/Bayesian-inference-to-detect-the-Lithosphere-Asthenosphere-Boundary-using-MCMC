function [E,Cov,sampled_m_IS] = momentum_by_ISMC(k,nSVI,nIS,priorInfo,my_likelihood,fwd)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generate samples following of a Markov chain whose stationary distribution
is the posterior distributiomn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    np: 
    nSVI:
    nIS: 
    priorPdf:
    my_likelihood:
Outputs: 
    E: Expectation         
    Cov: Covariance | dim: np x np x nIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}



% Stage 1: Find the optimal IS density (or trial density) from the parametric
% density family.
% Here we fix to the multivariate gaussian distribution
% generate samples acoording to the given distribution

% generate C samples acoording to the prior distribution
switch fwd.mtype
     case "1"
        sampled_B = priorInfo.bdist.random(nSVI,k-1);
        sampled_D = priorInfo.ddist.random(nSVI,k);
        % evaluate the likelihood , mlikelihood at the generated samples
        sum_mleval = 0;
        Leval_SVI = zeros(nSVI,1);
        for i  = 1:nSVI
             if rem(i,10) == 1
                fprintf("IS sample : %d \n",i);
             end
            leval_SVI = my_likelihood(sampled_B(i,:),sampled_D(i,:));
            Leval_SVI(i,1) = leval_SVI;
            sum_mleval = sum_mleval + [sampled_B(i,:)';sampled_D(i,:)']*leval_SVI; %  np x nSVI
        end
        mu_star = sum_mleval/ sum(Leval_SVI,1); % np x 1
        mmLeval_SVI = zeros(nSVI,2*k-1,2*k-1); % np x np x SVI
        for I = 1:nSVI
            mmLeval_SVI (I,:,:) = ([sampled_B(I,:)';sampled_D(I,:)']-mu_star)*([sampled_B(I,:)';sampled_D(I,:)']-mu_star)'* Leval_SVI(I,1);
        end

        covMatrix_star  = squeeze(sum(mmLeval_SVI,1)/sum(Leval_SVI,1)) + eps; % to make sure it is positive definite

        % Stage 2: generate samples according to the optimal IS density
        sampled_m_IS = mvnrnd(mu_star,covMatrix_star,nIS)'; % np x nIS

        % evaluate the likelihood, mLikelihood at the generated samples
        %sum_mleval_IS = 0;
        Leval_IS = zeros(nIS,1);
        mLeval_IS = zeros(nIS,2*k-1);
        for i = 1:nIS
            leval_IS = my_likelihood(sampled_m_IS(1:k-1,i),sampled_m_IS(k:end,i));
            Leval_IS(i,1) = leval_IS;
            mLeval_IS(i,:) = sampled_m_IS(:,i).*leval_IS;
        end
        % evaluate the prior density at the generated points
        Prior_eval = priorInfo.piPdf(sampled_m_IS(1:k-1,:),sampled_m_IS(k:end,:));% nIS x 1
        % evaluate the importance sampling weights the the generated samples
        IS_eval = mvnpdf(sampled_m_IS',mu_star',covMatrix_star);% nIS x 1
        Weights = Prior_eval./IS_eval;% nIS x 1
        E = mLeval_IS'*Weights /(Leval_IS'*Weights); % np x 1

        % weights , more samples nIS
        mmLeval_IS = zeros(nIS,2*k-1,2*k-1);
        for I = 1:nIS
            mmLeval_IS (I,:,:) = (sampled_m_IS(:,I)-E)*(sampled_m_IS(:,I)-E)'* Leval_IS(I,1)*Weights(I);
        end
        Cov = squeeze(sum(mmLeval_IS,1)/(Leval_IS'*Weights));
    case "2"
        sampled_C = priorInfo.cdist.random(nSVI,k);
        sampled_D = priorInfo.ddist.random(nSVI,k);
        %sampled_m_SVI = mvnrnd(priorInfo.mu,priorInfo.covMatrix,nSVI)'; % np x nSVI

        % evaluate the likelihood , mlikelihood at the generated samples
        sum_mleval = 0;
        Leval_SVI = zeros(nSVI,1);
        for i  = 1:nSVI
             if rem(i,10) == 1
                fprintf("IS sample : %d \n",i);
             end
            leval_SVI = my_likelihood(sampled_C(i,:),sampled_D(i,:));
            Leval_SVI(i,1) = leval_SVI;
            sum_mleval = sum_mleval + [sampled_C(i,:)';sampled_D(i,:)']*leval_SVI; %  np x nSVI
        end
        mu_star = sum_mleval/ sum(Leval_SVI,1); % np x 1
        mmLeval_SVI = zeros(nSVI,2*k,2*k); % np x np x SVI
        for I = 1:nSVI
            mmLeval_SVI (I,:,:) = ([sampled_C(I,:)';sampled_D(I,:)']-mu_star)*([sampled_C(I,:)';sampled_D(I,:)']-mu_star)'* Leval_SVI(I,1);
        end

        covMatrix_star  = squeeze(sum(mmLeval_SVI,1)/sum(Leval_SVI,1)) + eps; % to make sure it is positive definite
        %covMatrix_star2  = squeeze(sum(mmLeval_SVI2,1)/sum(Leval_SVI,1))-mu_star'*mu_star;

        % Stage 2: generate samples according to the optimal IS density
        sampled_m_IS = mvnrnd(mu_star,covMatrix_star,nIS)'; % np x nIS

        % evaluate the likelihood, mLikelihood at the generated samples
        %sum_mleval_IS = 0;
        Leval_IS = zeros(nIS,1);
        mLeval_IS = zeros(nIS,2*k);
        for i = 1:nIS
            leval_IS = my_likelihood(sampled_m_IS(1:k,i),sampled_m_IS(k+1:end,i));
            Leval_IS(i,1) = leval_IS;
            %sum_mLeval_IS = sum_mLeval_IS + sampled_m_IS(:,i).*leval_IS;
            mLeval_IS(i,:) = sampled_m_IS(:,i).*leval_IS;
        end
        % evaluate the prior density at the generated points
        Prior_eval = priorInfo.piPdf(sampled_m_IS(1:k,:),sampled_m_IS(k+1:end,:));% nIS x 1
        % evaluate the importance sampling weights the the generated samples
        %Prior_eval = mvnpdf(sampled_m_IS',priorInfo.mu',priorInfo.covMatrix);% nIS x 1
        IS_eval = mvnpdf(sampled_m_IS',mu_star',covMatrix_star);% nIS x 1
        Weights = Prior_eval./IS_eval;% nIS x 1
        E = mLeval_IS'*Weights /(Leval_IS'*Weights); % np x 1

        %E = sum(mLeval_IS,1) / sum(Leval_IS,1);

        % weights , more samples nIS
        mmLeval_IS = zeros(nIS,2*k,2*k);
        for I = 1:nIS
            mmLeval_IS (I,:,:) = (sampled_m_IS(:,I)-E)*(sampled_m_IS(:,I)-E)'* Leval_IS(I,1)*Weights(I);
        end
        Cov = squeeze(sum(mmLeval_IS,1)/(Leval_IS'*Weights));
    case "3"
        sampled_C = priorInfo.cdist.random(nSVI,k-1);
        sampled_D = priorInfo.ddist.random(nSVI,k);
        %sampled_m_SVI = mvnrnd(priorInfo.mu,priorInfo.covMatrix,nSVI)'; % np x nSVI

        % evaluate the likelihood , mlikelihood at the generated samples
        sum_mleval = 0;
        Leval_SVI = zeros(nSVI,1);
        for i  = 1:nSVI
             if rem(i,10) == 1
                fprintf("IS sample : %d \n",i);
             end
            leval_SVI = my_likelihood([fwd.c1;sampled_C(i,:)],sampled_D(i,:));
            Leval_SVI(i,1) = leval_SVI;
            sum_mleval = sum_mleval + [sampled_C(i,:)';sampled_D(i,:)']*leval_SVI; %  np x nSVI
        end
        mu_star = sum_mleval/ sum(Leval_SVI,1); % np x 1
        mmLeval_SVI = zeros(nSVI,2*k-1,2*k-1); % np x np x SVI
        for I = 1:nSVI
            mmLeval_SVI (I,:,:) = ([sampled_C(I,:)';sampled_D(I,:)']-mu_star)*([sampled_C(I,:)';sampled_D(I,:)']-mu_star)'* Leval_SVI(I,1);
        end

        covMatrix_star  = squeeze(sum(mmLeval_SVI,1)/sum(Leval_SVI,1)) + eps; % to make sure it is positive definite
        %covMatrix_star2  = squeeze(sum(mmLeval_SVI2,1)/sum(Leval_SVI,1))-mu_star'*mu_star;

        % Stage 2: generate samples according to the optimal IS density
        sampled_m_IS = mvnrnd(mu_star,covMatrix_star,nIS)'; % np x nIS

        % evaluate the likelihood, mLikelihood at the generated samples
        %sum_mleval_IS = 0;
        Leval_IS = zeros(nIS,1);
        mLeval_IS = zeros(nIS,2*k-1);
        for i = 1:nIS
            leval_IS = my_likelihood([fwd.c1;sampled_m_IS(1:k-1,i)],sampled_m_IS(k:end,i));
            Leval_IS(i,1) = leval_IS;
            %sum_mLeval_IS = sum_mLeval_IS + sampled_m_IS(:,i).*leval_IS;
            mLeval_IS(i,:) = sampled_m_IS(:,i).*leval_IS;
        end
        % evaluate the prior density at the generated points
        Prior_eval = priorInfo.piPdf(sampled_m_IS(1:k-1,:),sampled_m_IS(k:end,:));% nIS x 1
        % evaluate the importance sampling weights the the generated samples
        %Prior_eval = mvnpdf(sampled_m_IS',priorInfo.mu',priorInfo.covMatrix);% nIS x 1
        IS_eval = mvnpdf(sampled_m_IS',mu_star',covMatrix_star);% nIS x 1
        Weights = Prior_eval./IS_eval;% nIS x 1
        E = mLeval_IS'*Weights /(Leval_IS'*Weights); % np x 1

        % weights , more samples nIS
        mmLeval_IS = zeros(nIS,2*k-1,2*k-1);
        for I = 1:nIS
            mmLeval_IS (I,:,:) = (sampled_m_IS(:,I)-E)*(sampled_m_IS(:,I)-E)'* Leval_IS(I,1)*Weights(I);
        end
        Cov = squeeze(sum(mmLeval_IS,1)/(Leval_IS'*Weights));
end
