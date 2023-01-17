function er = Error(fwd,m)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Error_ob random variable definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    fwd: 
    m:              sample to be evaluated | dim: np x 1

Outputs: 
    er:                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
% if isempty(fwd.table)
%    sol = forward(fwd,c,d);
% else
%    sol = getSolutionFromTable(fwd,c,d);
% end
labmesh = setParameterization(m,fwd);
sol = forward(fwd,labmesh);
%sol2 = forward(fwd, m);%obs2 = fwd.B0 * sol2;
obs = fwd.B0*sol;
er = (obs - fwd.observedData);
end
