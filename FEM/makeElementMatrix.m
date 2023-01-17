function [Ke,Ge,fe] = makeElementMatrix(Xe,fwd,labmesh)
%
mesh  = fwd.mesh;
material = fwd.material;

[nOfVelNodesPerElement] = size(mesh.T,2);
[~,nOfPreNodesPerElement] = size(mesh.TP);
nOfVelDofPerElement = mesh.nsd * nOfVelNodesPerElement;
nOfPreDofPerElement = nOfPreNodesPerElement;

% initialisation
Ke = zeros(nOfVelDofPerElement,nOfVelDofPerElement);
Ge = zeros(nOfPreDofPerElement,nOfVelDofPerElement);
fe = zeros(nOfVelDofPerElement,1);

N = mesh.refElem.N;
Nxi = mesh.refElem.Nxi;
Neta = mesh.refElem.Neta;
NP = mesh.refElem.NP;
pespg = mesh.refElem.pespg;

% number of gauss points
ngaus = mesh.refElem.nGaussPoints;

% Get the global coordinates of gauss points based on the global coord of
% the velocity nodes
Xgaus =  N*Xe;
%LABmesh = setParameterization(b,d,fwd);
labDepth   = LABdepth(Xgaus(:,1),labmesh);
[g,rho,nu] = getProperties(material,labDepth,Xgaus(:,2));

for igaus = 1: ngaus
   jacob = [Nxi(igaus,:) *Xe(:,1)  Nxi(igaus,:) *Xe(:,2); ...
            Neta(igaus,:)*Xe(:,1)  Neta(igaus,:)*Xe(:,2)];
   dvolu = pespg(igaus)*det(jacob);
   % shape functions derivatives in cartesian coordinates
   res = jacob\[Nxi(igaus,:);Neta(igaus,:)];
   % gradient of the shape functions (cartesian coords)
   Nx = [reshape([1;0]*res(1,:),1,nOfVelDofPerElement); ...
         reshape([0;1]*res(1,:),1,nOfVelDofPerElement)];
   Ny = [reshape([1;0]*res(2,:),1,nOfVelDofPerElement); ...
         reshape([0;1]*res(2,:),1,nOfVelDofPerElement)];
   % divergence of shape functions
   dN = reshape(res,1,nOfVelDofPerElement);
   % element matrices
   Ke = Ke + nu(igaus)  * (Nx'*Nx + Ny'*Ny) * dvolu;
   Ge = Ge - NP(igaus,:)' * dN * dvolu;
   N_igaus = reshape([0;1]*N(igaus,:),1,nOfVelDofPerElement);
   fe = fe + N_igaus'*(rho(igaus) *g*dvolu);
end
end


