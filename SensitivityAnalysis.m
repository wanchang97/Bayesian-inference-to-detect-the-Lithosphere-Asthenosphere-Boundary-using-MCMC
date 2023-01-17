% sensitivity analysis
% we will calculate (v(c+dc)-v(c))/dc and (v(d+dd)-v(d))/dd

%% The basic settings
clearvars
close all
clc
%%
addpath("mesh");
addpath("utils");
addpath("inverse");
addpath("FE");
addpath("plot");
figurepath = 'figuresVoronoi';
%% Forward problem
% space dimension
fwd.nsd = 2;
% Domain
fwd.geometry.lx = 400*10^3; % m
fwd.geometry.lz = 660*10^3; % m

% Physical properties
fwd.material.gravity = 9.81;
fwd.material.nu  = [10^23,10^20];% Pa∙s; viscocity in omega1 and omega2
fwd.material.rho = [3300, 3280]; % kg/m^3; density in omega1 and omega2nu

% Discretization
% element type: 1 = squares, 2 = triangles
FE.elemType = 1;
% number of nodes per element
FE.nenVel = 9; FE.nenPre = 4;
% number of elements in each direction
FE.nx = 30; FE.nz = 30;
FE.nGaussPoints = 9;
fwd.mesh = createMeshForVelocityandPressure(fwd.geometry,FE);

%% Observables
sensorExactLocations = ...
    [fwd.geometry.lx/3,   fwd.geometry.lz*2/5  ; ...
    fwd.geometry.lx*1/4, fwd.geometry.lz*2/5  ; ...
    fwd.geometry.lx*1/5, fwd.geometry.lz*9/20; ...
    fwd.geometry.lx/4,   fwd.geometry.lz*9/20  ; ...
    fwd.geometry.lx/2,   fwd.geometry.lz*9/20  ];
nSensorLocation = size(sensorExactLocations,1);
% Sensors (at the closest node)
for I = 1:size(sensorExactLocations,1)
    sensor(2*I-1) = setSensor(fwd.mesh,sensorExactLocations(I,:),1); % x
    sensor(2*I)   = setSensor(fwd.mesh,sensorExactLocations(I,:),2); % y
end

% Quantity of interest
% qoi_scaling_constant = 1E11;
fwd.B0 = makeQoIMatrix(sensor,fwd.mesh);

%% True parameters (to be recovered)
% Parameterization fixe one input parameter
k = 2; true_d = fwd.geometry.lz*[1/3;1/2];true_c = fwd.geometry.lx*[1/4;1/2];
%k = 3; true_d = fwd.geometry.lz*[1/3;1/2;2/5];true_c = fwd.geometry.lx*[1/4;1/2;2/3];
fwd.k = k;fwd.mtype = "3";
switch fwd.mtype
    case "2"
        fwd.np = fwd.k + fwd.k * (fwd.nsd-1);
    case "3"
        fwd.np = fwd.k + (fwd.k-1) * (fwd.nsd-1);
end

%% the true value
true_m = [true_c; true_d];
fwd.c1 = true_c(1);% only for 1D case
true_LABmesh = setParameterization(true_c,true_d,fwd);
sol = forward(fwd,true_c,true_d);
v = sol(1:fwd.mesh.nOfVelDof);
v = reshape(v',2,[])';
figure,clf
quiver(fwd.mesh.X(:,1),fwd.mesh.X(:,2),v(:,1),v(:,2));
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
axis equal tight
plotLAB(true_LABmesh,true_d)
hold on
for i = 1:nSensorLocation
    plot(sensorExactLocations(i,1),sensorExactLocations(i,2),'g*');
end
%% Observations
true_observation = fwd.B0*sol;
fwd.observedData = true_observation; % no noise at the moment


%% The alteration 1 when change the coordinates of the sites
dc2 = 0.5 * true_c(2)
new_c2 = [fwd.c1;true_c(2) + dc2];
new_m_c2 = [new_c2;true_d];
new_LABmesh_c2 = setParameterization(new_c2,true_d,fwd);
new_sol_c2 = forward(fwd,new_LABmesh_c2.C,new_LABmesh_c2.d);
new_observation_c2 = fwd.B0*new_sol_c2;
sensitivityC2 = norm(new_observation_c2-true_observation)/dc2
norm(new_sol_c2-sol);
%% 
dd1 = 0.5 * true_d(1)
new_d1 = [true_d(1) + dd1;true_d(2)];
new_m_d1 = [true_c;new_d1];
new_LABmesh_d1 = setParameterization(true_c,new_d1,fwd);
new_sol_d1 = forward(fwd,new_LABmesh_d1.C,new_LABmesh_d1.d);
new_observation_d1 = fwd.B0*new_sol_d1;
sensitivity_d1 = norm(new_observation_d1-true_observation)/dd1
norm(new_sol_d1-sol);
%%
dd2 = 0.5 * true_d(2)
new_d2 = [true_d(1) ;true_d(2)+dd2];
new_m_d2 = [true_c;new_d2];
new_LABmesh_d2 = setParameterization(true_c,new_d2,fwd);
new_sol_d2 = forward(fwd,new_LABmesh_d2.C,new_LABmesh_d2.d);
new_observation_d2 = fwd.B0*new_sol_d2;
sensitivity_d2 = norm(new_observation_d2-true_observation)/dd2
norm(new_sol_d2-sol);