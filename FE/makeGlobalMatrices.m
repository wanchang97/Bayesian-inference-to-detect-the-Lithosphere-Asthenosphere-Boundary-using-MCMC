function [K,G,f] = makeGlobalMatrices(fwd,labmesh)
%
mesh  = fwd.mesh;

X  = mesh.X;
T  = mesh.T;
TP = mesh.TP;

% number of nodes per element different in velocity and pressure
[nOfElements,nOfVelocityNodesPerElement] = size(T);
[~,nOfPressureNodesPerElement] = size(TP);

% number of dof per element
nOfVelocityDofPerElement = mesh.nsd * nOfVelocityNodesPerElement;
nOfPressureDofPerElement = nOfPressureNodesPerElement;

% element matrix size
mKe = nOfVelocityDofPerElement;
nKe = nOfVelocityDofPerElement;
mGe = nOfPressureDofPerElement;
nGe = nOfVelocityDofPerElement;
% initialisation
allK = zeros(mKe*nKe,nOfElements);
allKI = zeros(mKe*nKe,nOfElements);
allKJ = zeros(mKe*nKe,nOfElements);
allG = zeros(mGe*nGe,nOfElements);
allGI = zeros(mGe*nGe,nOfElements);
allGJ = zeros(mGe*nGe,nOfElements);
f = zeros(mesh.nOfVelDof,1);

% loop in elements
for ielem = 1:nOfElements
    % velocity nodes ux, uy 
    Te = reshape([2*T(ielem,:)-1;2*T(ielem,:)],1,nOfVelocityDofPerElement);
    % pressure nodes
    TeP = TP(ielem,:);
    % global coords of the velocity nodes
    Xe = X(T(ielem,:),:);
    % element matrix
    [Ke,Ge,fe] = makeElementMatrix(Xe,fwd,labmesh);

    % store elemental matrices
    allK(:,ielem) = Ke(:);
    [mi,mj] = meshgrid(Te,Te);
    allKI(:,ielem) = mi(:);
    allKJ(:,ielem) = mj(:);

    allG(:,ielem) = Ge(:);
    [mi,mj] = meshgrid(Te,TeP);
    allGI(:,ielem) = mi(:);
    allGJ(:,ielem) = mj(:);
    % assemby force term
    f(Te,1) = f(Te,1) + fe;
end

K = sparse(allKJ,allKI,allK);
G = sparse(allGJ,allGI,allG);
end

