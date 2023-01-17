function [C,D] = RW_Bodin(piPdf,likelihood,c0,d0,nRW,sigmaC,sigmad)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generate samples following of a Markov chain whose stationary distribution
is the posterior distribution
This will distinguish c and d and update them one by one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    piPdf           : 
    likelihood:     :
    c0              :The initial sites location vector |dim n x nsd
    d0              :The initial LAB depth vector |dim n x 1
    nRW             :The number of maximal interations
    sigmaC          :sqrt of diagonal of the covariance matrix|dim :n x 1
    sigmad          :sqrt of 
Outputs: 
    C               :Samples of sites in Markov chain |dim: nsamples x (n x nsd)
    D               :Samples of depth in Markov chain |dim: nsamples x (n x 1)
    M               :Samples of the markov chain |dim: nsamples x (n x (nsd+1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%}
nsd = size(c0,2);
n = size(c0,1);
C = zeros(nRW,n,nsd);
D = zeros(nRW,n,1);
c_candidate = c0;
d_candidate = d0;
m_candidate = [c0,d0];
C(1,:,:) = c0;
D(1,:,:) = d0;
muCRW = zeros(nsd,1);
covCRW = diag(sigmaC).^2;
muDRW = zeros(1,1);
covDRW = diag(sigmad).^2;
for i = 2:nRW
    % pick up a cell randomly
    cell = randi([1,n],1,1);
    % either update its position 
    if rem(i-1,2) % for even iteration
        % propose a new candidate coordinate from the proposal density
        %z_coordinate = normrnd(muCRW,covCRW,1,nsd);
        z_coordinate = mvnrnd(muCRW,covCRW);
        % How to make the new candidate inside the domain?
        %z_coordinate = mvnrnd(muCRW,covCRW,1,nsd);
        c_candidate_cell = C(i-1,cell,:) + z_coordinate;
        c_candidate(cell,:) = c_candidate_cell;
    else
    % or update its depth
        z_depth = mvnrnd(muDRW,covDRW);
        %sz_depth = mvnrnd(muDRW,covDRW,1,1);
        d_candidate_cell = D(i-1,cell,:) + z_depth;
        d_candidate(cell,:) = d_candidate_cell;
    end
    % solve the forward problem and the convergence ratio
    %r = prior(c_candidate,d_candidate)*likelihood(c_candidate,d_candidate)/(prior(C(i-1,:,:),D(i-1,:,:))*likelihood(C(i-1,:,:),D(i-1,:,:)));
    % for the uniform prior distribution 
    r = likelihood(c_candidate,d_candidate)/likelihood(C(i-1,:,:),D(i-1,:,:));
    alpha = min(1,r);
    % generate a sample u from the uniform distribution
    u = unifrnd(0,1,1);
    if u <= alpha % accept the sample
        C(i,:,:) = c_candidate;
        D(i,:,:) = d_candidate;
    else % reject the sample
        C(i,:,:) = C(i-1,:,:);
        D(i,:,:) = D(i-1,:,:);
    end
end
% the burn-in period is set to be 0.1 *nRW
nburn = nRW * 0.1;
C = C(nburn:nRW,:,:);
D = D(nburn:nRW,:,:);
end