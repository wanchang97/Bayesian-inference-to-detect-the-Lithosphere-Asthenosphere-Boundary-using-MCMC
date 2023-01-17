function [B,D,acc] = RW(piPdf,likelihood,b0,d0,nRW,sigmab,sigmad)
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
    sigmab          :sqrt of diagonal of the covariance matrix|dim :n x 1
    sigmad          :sqrt of 
Outputs: 
    B               :Samples of sites in Markov chain |dim: nsamples x (n x nsd)
    D               :Samples of depth in Markov chain |dim: nsamples x (n x 1)
    acc             :The number of samples being accepted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
nsd = size(b0,2);
k = size(d0,1);

B = zeros(nRW,k-1,nsd);
D = zeros(nRW,k,1);
B(1,:) = b0;
D(1,:) = d0;
muBRW = zeros((k-1)*nsd,1);
covBRW = diag(sigmab.^2*ones(1,k-1));
muDRW = zeros(k,1);
covDRW = diag(sigmad.^2*ones(1,k));
acc = 0;
for i = 2:nRW
    if rem(i,10) == 1
        fprintf('RW sample: %d \n ', i )
    end
    % propose a new candidate coordinate from the proposal density
    %z_coordinate = normrnd(muCRW,covCRW,1,nsd);
    z_coordinate = mvnrnd(muBRW,covBRW)';
    % How to make the new candidate inside the domain?
    %z_coordinate = mvnrnd(muCRW,covCRW,1,nsd);
    b_candidate = B(i-1,:)' + z_coordinate;
    % and update its depth
    z_depth = mvnrnd(muDRW,covDRW)';
    %sz_depth = mvnrnd(muDRW,covDRW,1,1);
    d_candidate = D(i-1,:)' + z_depth;
    % solve the forward problem and the convergence ratio
    %r = prior(c_candidate,d_candidate)*likelihood(c_candidate,d_candidate)/(prior(C(i-1,:,:),D(i-1,:,:))*likelihood(C(i-1,:,:),D(i-1,:,:)));
    % for the uniform prior distribution
    r = likelihood(b_candidate,d_candidate)*piPdf(b_candidate,d_candidate)/(likelihood(B(i-1,:)',D(i-1,:)')*piPdf(B(i-1,:)',D(i-1,:)'));
    alpha = min(1,r);
    % generate a sample u from the uniform distribution
    u = unifrnd(0,1,1);
    if u <= alpha % accept the sample
        B(i,:) = b_candidate;
        D(i,:) = d_candidate;
        acc = acc + 1;
    else % reject the sample
        B(i,:) = B(i-1,:);
        D(i,:) = D(i-1,:);
    end
end
% the burn-in period is set to be 0.1 *nRW
%nburn = nRW * 0.1;
%C = C(nburn:nRW,:,:);
%D = D(nburn:nRW,:,:);
end