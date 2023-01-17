function [E,Cov] = momentum_by_integrals(prior,my_likelihood,grid,np)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Evaluate the E and V using the Middle Points rule. however, with the limit
of our model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    prior:           pi = @(m) priorPdf(m,prior.mu,prior.covMatrix, prior.type);
    my_likelihood:   my_likelihood = @(m) likelihood(fwd,m,Covmatrix_obs);
    grid: 
    np:
Outputs: 
    E:                Expectation of Samples of the markov chain | dim: np x 1
    Cov:              Covariance matrix| dim: np x np 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
% definition of the quadrature
w = grid(2:end) - grid(1:end-1);
z = (grid(2:end) + grid(1:end-1)) / 2;
nOfIntegrationPoints = length(w);% 61

C = 0;
N = 0;
n = nOfIntegrationPoints^np;
s = nOfIntegrationPoints * ones(1,np);

storedLikelihoods = zeros(n,1);
ix = cell(np,1);
for I = 1:n
    [ix{:}] = ind2sub(s,I);
    m = z([ix{:}]);
    weight = prod(w([ix{:}]));
    %temp =  my_likelihood(m)
    storedLikelihoods(I) = my_likelihood(m);
    aux = storedLikelihoods(I)*prior(m)*weight;
    C = C + aux;
    N = N + m*aux;
end
E = N/C;

N2 = 0;
for I = 1:n
    [ix{:}] = ind2sub(s,I);
    m = z([ix{:}]);
    weight = prod(w([ix{:}]));
    aux = storedLikelihoods(I)*prior(m)*weight;
    N2 = N2 + (m-E)*(m-E)'*aux;
end
Cov = N2/C;

% switch np
%     case 2
%         % FIXME: only work for two parameters
%         storedLikelihoods = zeros(nOfIntegrationPoints,nOfIntegrationPoints);
%         for I1 = 1:nOfIntegrationPoints
%             for I2 = 1:nOfIntegrationPoints
%                 m = [z(I1); z(I2)];
% 
%                 storedLikelihoods(I1,I2) = my_likelihood(m);
%                 aux = storedLikelihoods(I1,I2) * prior(m) * w(I1) * w(I2);
% 
%                 C = C + aux;
%                 N = N + m*aux;
%             end
%         end
%         E = N / C;
% 
%         % FIXME: 2 parameter case
%         N2 = 0;
%         for I1 = 1:nOfIntegrationPoints
%             for I2 = 1:nOfIntegrationPoints
%                 m = [z(I1); z(I2)];
%                 aux = storedLikelihoods(I1,I2) * prior(m) * w(I1) * w(I2);
%                 N2 = N2 + (m-E)*(m-E)'*aux;
%             end
%         end
%         Cov = N2 / C;
%     case 3
% % FIXME: only work for 3 parameters
% storedLikelihoods = zeros(nOfIntegrationPoints,nOfIntegrationPoints,nOfIntegrationPoints);
%         for I1 = 1:nOfIntegrationPoints
%             for I2 = 1:nOfIntegrationPoints
%                 for I3 = 1:nOfIntegrationPoints
%                 m = [z(I1); z(I2); z(I3)];
% 
%                 storedLikelihoods(I1,I2,I3) = my_likelihood(m);
%                 aux = storedLikelihoods(I1,I2,I3) * prior(m) * w(I1) * w(I2) * w(I3);
% 
%                 C = C + aux;
%                 N = N + m*aux;
%                 end
%             end
%         end
%         E = N / C;
%         % FIXME: 3 parameter case
%         N2 = 0;
%         for I1 = 1:nOfIntegrationPoints
%             for I2 = 1:nOfIntegrationPoints
%                 for I3 = 1:nOfIntegrationPoints
%                 m = [z(I1); z(I2); z(I3)];
%                 aux = storedLikelihoods(I1,I2,I3) * prior(m) * w(I1) * w(I2) *w(I3);
%                 N2 = N2 + (m-E)*(m-E)'*aux;
%                 end
%             end
%         end
%         Cov = N2 / C;
% end
% 
% 
