function [E,Cov,CIdelta] = Ensemble_RW(M,fwd,nburn,t,alpha)
%{
Input:
    B:              samples of depth
    D:              samples of depth
Output: 
    E:              E = int x f(x) dx approx E_hat = mean (x_i)
    Cov:            Cov = int (x-E)(x-E)^T f(x) dx approx Cov_hat = mean
    (x_i-E)(x_i-E)^T
E_uncertainty:      

%}
k = fwd.k;
nsamples = size(M,1)-nburn;
switch fwd.mtype
    case 0
        M = M(nburn:t:end,:,:);
        E = mean(M,1);
        Cov = transpose(M-E)*(M-E)/nsamples;
        CIdelta = CIinterval(alpha,M);
        E = E';

    case 1
        M = M(nburn:t:end,:,:);
        B = M(:,1:k-1);
        D = M(:,k:end);
        EB = mean(B,1);
        ED = mean(D,1);
        E = [EB,ED];
        Cov = transpose(M-E)*(M-E)/nsamples;
        CIdelta = CIinterval(alpha,M);
        E = E';

end
end
