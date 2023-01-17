function [Mdir,Vdir] = bcCavity(X,nVelocityDof)
%{
---------------------------------------------------------------------------
Builds the matrix Adir, and vector Vdir to impose the boundary conditions
of the cavity problem using the lagrange multipliers.
Assume a square domain
---------------------------------------------------------------------------
Input: 
    X         : The coordinate of the velocity nodes
nVelocityDof  : number of velocity DDF per element
      
Output: 
    Mdir      : element stiffness matrix K dimension nVelocityDofPerElement x nVelocityDofPerElement 
    Vdir      : element gradient matrix Gs with dimension nPressureNodesPerElement x nVelocityDofPerElement 
---------------------------------------------------------------------------
%}
% find domain size
xmin = min(X(:,1));
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax= max(X(:,2));
% find the nodes on each boundary using tolerance
tol = 1e-6;
% boundary nodes except the corner
topNodes = reshape(find(abs(X(:,2) - ymin) <= tol & xmin < X(:,1) & X(:,1) < xmax),[],1);
botNodes = reshape(find(abs(X(:,2) - ymax) <= tol & xmin < X(:,1) & X(:,1) < xmax),[],1);
lefNodes = reshape(find(abs(X(:,1) - xmin) <= tol & ymin < X(:,2) & X(:,2) < ymax),[],1);
rigNodes = reshape(find(abs(X(:,1) - xmax) <= tol & ymin < X(:,2) & X(:,2) < ymax),[],1);
% Corner nodes
lefuppercorNode       = find(abs(X(:,1)-xmin) <= tol & abs(X(:,2)-ymin) <= tol);
lefbottomcorNode      = find(abs(X(:,1)-xmin) <= tol & abs(X(:,2)-ymax) <= tol);
rightuppercorNode     = find(abs(X(:,1)-xmax) <= tol & abs(X(:,2)-ymin) <= tol);
rightbottomcorNode    = find(abs(X(:,1)-xmax) <= tol & abs(X(:,2)-ymax) <= tol);
corNodes = reshape([lefuppercorNode,lefbottomcorNode,rightuppercorNode,rightbottomcorNode],[],1);
% build matrix: [node  bc_value_on_that_node] normal velocity is 0
nn = length(topNodes);
bctop = [2*topNodes zeros(nn,1)];% only need to assign the bottom part to be zero and the above part arbitrary

nn = length(botNodes);
bcbot = [2*botNodes zeros(nn,1)];

nn = length(lefNodes);
bclef = [2*lefNodes-1 zeros(nn,1)];

nn = length(rigNodes);
bcrig = [2*rigNodes-1 zeros(nn,1)];

nn = length(corNodes);
bccor = [[2*corNodes-1;2*corNodes],[zeros(2*nn,1)]];
% all boundary conditions together
bc = [bctop; bcbot; bclef; bcrig;bccor];

% matrix to impose BC using lagrange multipliers
nDirichletBC = size(bc,1);
Mdir = zeros(nDirichletBC,nVelocityDof);
Mdir(:,bc(:,1)) = eye(nDirichletBC);
% vector to impose BC using lagrange multipliers
Vdir = bc(:,2);
end

