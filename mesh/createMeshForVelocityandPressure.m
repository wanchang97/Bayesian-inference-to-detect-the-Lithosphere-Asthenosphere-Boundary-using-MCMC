function theMesh = createMeshForVelocityandPressure(geometry,FE)
%

% Element type: 1= squares, 2= triangles
elemType = FE.elemType;
% number of elements in each direction
nx = 2*FE.nx; nz = 2*FE.nz;
% number of nodes per element
nenVel = FE.nenVel; nenPre = FE.nenPre;
% domain size
xmin = 0; xmax = geometry.lx;
zmin = 0; zmax = geometry.lz;
% create mesh for velocity
[theMesh.X,theMesh.T] = createVelocityMesh(elemType,nenVel,xmin,xmax,zmin,zmax,nx,nz);
% create mesh for pressure
[theMesh.XP,theMesh.TP] = createPressureMesh(elemType,nenPre,theMesh.X,theMesh.T,nx,nz);
theMesh.elemType = elemType;
theMesh.nsd = 2; % number of space dimensions

theMesh.nx = FE.nx; 
theMesh.nz = FE.nz;
theMesh.nOfVelDof = size(theMesh.X,1) * theMesh.nsd;
theMesh.nOfPreDof = size(theMesh.XP,1);


% Reference Element
refElem.nGaussPoints = FE.nGaussPoints;
[refElem.pospg,refElem.pespg] = quadrature(elemType,refElem.nGaussPoints);
[refElem.N,refElem.Nxi,refElem.Neta] = ...
   shapeFunctions(elemType,FE.nenVel,refElem.pospg);
[refElem.NP,~,~] = shapeFunctions(elemType,FE.nenPre,refElem.pospg);

theMesh.refElem = refElem;
end