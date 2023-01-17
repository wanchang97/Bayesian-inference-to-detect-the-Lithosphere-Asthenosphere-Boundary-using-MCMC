function logL = loglikelihood(fwd,m,covmatrix)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logLikelihood definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    fwd : 
    b   :         width of the column | dim: (k-1)x 1
    d   :         depth of the column | dim: k x 1
    covmatrix: 
Outputs: 
    logL:                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

n_ob = size(fwd.B0,1);
%k = fwd.k;
logcte = -1/2 * (n_ob*log(2*pi)+sum(log(diag(covmatrix))));
switch fwd.mtype
    case 0
        dmin = fwd.param.dmin;
        dmax = fwd.param.dmax;
        d = m;
        outOfRange = any(d < dmin)| any(d > dmax);
    case 1
        bmin = fwd.param.bmin;
        bmax = fwd.param.bmax;
        dmin = fwd.param.dmin;
        dmax = fwd.param.dmax;
        %b = m(1:k-1);
        %d = m(k:end);
        LABmesh = setParameterization(m,fwd);
        BoundaryLoc = LABmesh.X(2:end);
        d = LABmesh.d;
        outOfRange = any(BoundaryLoc < bmin) | any(BoundaryLoc > bmax)| any(d < dmin)| any(d > dmax);
end
if outOfRange
    logL = -10^7; % -Inf or 0?

else
    err  = Error(fwd,m);
    logL = -1/2 * err'* (covmatrix\err)+logcte;
end
end


