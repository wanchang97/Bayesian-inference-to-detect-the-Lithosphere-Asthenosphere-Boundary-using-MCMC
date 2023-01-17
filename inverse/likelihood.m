function L = likelihood(fwd,m,covmatrix)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Likelihood definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    fwd: 
   % C:               samples to be evaluated | dim: 2np x nsamples
   % D:               samples to be evaluated | dim: np x nsamples
    covmatrix: 
Outputs: 
    L:                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
k = fwd.k;
n_ob = size(fwd.B0,1);
switch fwd.mtype
    case 0
        dmin = fwd.param.dmin;
        dmax = fwd.param.dmax;%
        cte = 1/(sqrt(det(covmatrix)*(2*pi).^n_ob));
        d = m;
        outOfRange =  any(d < dmin)| any(d > dmax);
    case 1
        bmin = fwd.param.bmin;
        bmax = fwd.param.bmax;
        dmin = fwd.param.dmin;
        dmax = fwd.param.dmax;%
        cte = 1/(sqrt(det(covmatrix)*(2*pi).^n_ob));
        b = m(1:k-1);
        d = m(k:end);
        outOfRange = any(b < bmin) | any(b > bmax)| any(d < dmin)| any(d > dmax);
end
if outOfRange
    L= 0;
else
    err  = Error(fwd,m);
    L = cte * exp(-1/2 * err'* (covmatrix\err));
end
end


