function sol = forward(fwd,labmesh)
% m dim np x 1
%% Build global matrix
[K,G,f] = makeGlobalMatrices(fwd,labmesh);

%% Boundary conditions and make the solution unique by fixing one row
ksize = size(K,1); % ksize is the number of velocity dof
[Adir,bdir] = bcCavity(fwd.mesh.X,ksize);
lsize = size(Adir,1);% lsize is the number of boundary conditions, we need at least three in 2D and 6 in 3D to reduce the rank defeciency of K
G(1,:) = []; % To delete the first row of G equivalent to imposing p1 = 0
% build the global system
gsize = size(G,1); % now gsize+ 1 is the number of pressure dof
Z1 = zeros(lsize,lsize);
Z2 = zeros(gsize,lsize);
Z3 = zeros(gsize,gsize);
Atot  =[K    Adir' G'; ...
        Adir Z1    Z2'; ...
        G    Z2    Z3];
btot = [f; bdir;zeros(gsize,1)];

%% Solve the problem
aux = Atot\btot;

%% Extract the solutions
% extract the velocity
velocity  = aux(1:ksize);

% extract the pressure
ixs = ksize + lsize + 1;
pressure = [0;aux(ixs:end)];
sol = [velocity; pressure];


