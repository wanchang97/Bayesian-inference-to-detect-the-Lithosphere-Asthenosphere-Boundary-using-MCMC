function  C_emp = empiricalCovariance(Xt)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Evaluate an empirical covariance matrix based on available sample points
generated up to step t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    Xt        : all samples generated until step t | dim: t x (k-1) or t x k
Outputs: 
    C_emp     : empirical covarance matrix of time  t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
t = size(Xt,1);
k = size(Xt,2);

mux = mean(Xt,1);
covx = zeros(k,k); 

% Probably we dont need to calculate every time from the beginning?
for I = 1: t
    temp = transpose(Xt(I,:)-mux)*(Xt(I,:)-mux);
    covx = covx + temp;
end
C_emp = covx/(t-1);

C_emp_ref = cov(Xt,Xt);

end