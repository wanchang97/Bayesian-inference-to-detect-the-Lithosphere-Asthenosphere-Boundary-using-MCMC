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
%k = 2; true_d = fwd.geometry.lz*[1/3;7/12];true_c = fwd.geometry.lx*[1/4;1/2];
k = 3; true_d = fwd.geometry.lz*[1/3;1/2;2/5];true_c = fwd.geometry.lx*[1/4;1/2;2/3];
fwd.k = k;fwd.mtype = "2";
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
fwd.true_LABmesh = true_LABmesh;
sol = forward(fwd,true_LABmesh.C,true_LABmesh.d);
v = sol(1:fwd.mesh.nOfVelDof);
v = reshape(v',2,[])';
figure,clf
quiver(fwd.mesh.X(:,1),fwd.mesh.X(:,2),v(:,1),v(:,2));
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
axis equal tight
plotLAB(true_LABmesh,'plotSites',1)
hold on
for i = 1:nSensorLocation
    plot(sensorExactLocations(i,1),sensorExactLocations(i,2),'g*');
end
fileName = ['True LAB plot and np = ' num2str(fwd.np) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);

%% Observations
true_observation = fwd.B0*sol;
fwd.observedData = true_observation; % no noise at the moment

%% Backward problem

% parametric domains
fwd.param.dmin = 0.1 * fwd.geometry.lz;
fwd.param.dmax = 0.9 * fwd.geometry.lz;
fwd.param.cmin = 0.1 * fwd.geometry.lx;
fwd.param.cmax = 0.9 * fwd.geometry.lx;


DataName = ['fwd_' num2str(k) 'k_type' num2str(fwd.mtype)];
save(DataName,'fwd')
% Prior
priorInfo.ctype = 'uniform';
c_dist = Dist(priorInfo.ctype ,'PAR',[fwd.param.cmin,fwd.param.cmax]);
pic = @(c) c_dist.pdf(c);

priorInfo.dtype = 'uniform';
d_dist = Dist(priorInfo.dtype ,'PAR',[fwd.param.dmin,fwd.param.dmax]);
pid = @(d) d_dist.pdf(d);
mean_prior_d = d_dist.mean*ones(k,1);
priorInfo.cdist = c_dist;
priorInfo.ddist = d_dist;

switch fwd.mtype
    case "2"
        %fwd.np = (fwd.nsd-1)*fwd.k + fwd.k;
        mean_prior_c = c_dist.mean*ones(k,1);
        piPdf = @(c,d)priorPdf([c;d],[c_dist,d_dist]);
        priorInfo.piPdf = piPdf;
        LABmesh_prior = setParameterization(mean_prior_c,mean_prior_d,fwd);
        sol_prior = forward(fwd,mean_prior_c,mean_prior_d);
        v_prior = sol_prior(1:fwd.mesh.nOfVelDof);
        v_prior = reshape(v_prior',2,[])';
    case "3"
        %fwd.np = (fwd.nsd-1)*(fwd.k-1) + fwd.k; % totally 2*fwd.np for 2D model
        mean_prior_c = c_dist.mean*ones(k-1,1);
        piPdf = @(c,d)priorPdf([c;d],[c_dist,d_dist]);
        priorInfo.piPdf = piPdf;
        LABmesh_prior = setParameterization([fwd.c1;mean_prior_c],mean_prior_d,fwd);
        sol_prior = forward(fwd,[fwd.c1;mean_prior_c],mean_prior_d);
        v_prior = sol_prior(1:fwd.mesh.nOfVelDof);
        v_prior = reshape(v_prior',2,[])';
end
% Likelihood
sigma_ob_ratio = 0.2;
sigma_ob = abs(true_observation) * sigma_ob_ratio;
Covmatrix_obs = diag(sigma_ob.^2);
my_likelihood_fixsigma = @(c,d) likelihood(fwd,c,d,Covmatrix_obs);
%% load the data
fwddataName = ['fwd_' num2str(k) 'k_type' num2str(fwd.mtype)];
load(fwddataName);
LABdataName = ['LAB_RW_Result_' num2str(k) 'k_type' num2str(fwd.mtype)];
load(LABdataName);

%% Calcualte the stable chains
nburn = 0.2*nRW;
t = 5;
C_stable = C(nburn:t:end,:);
D_stable = D(nburn:t:end,:);

%% Plot the LAB tomography for each sample
f = figure
for i = 1:5
    LABmesh = setParameterization(C_stable(i,:),D_stable(i,:),fwd);
    subplot(1,6,i)
    plotLAB(LABmesh,'plotSites',1);
    hold on
    plotMesh(fwd.mesh.X,fwd.mesh.T,'elementType',1);
    title('Sample ',i)
    xlabel('x');ylabel('z');
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis equal tight
end
subplot(1,6,6)
solutionMap(C_stable,D_stable,fwd)
title('Average')
x0=10;
y0=10;
width=1000;
height=350;
set(gcf,'position',[x0,y0,width,height])
fileName = ['LAB plot for first five samples with average solution nRW = ' num2str(nRW) ' and np = ' num2str(fwd.np) ' and acc = ' num2str(acc) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
%% Plot the average solution map
solutionMap(C_stable,D_stable,fwd)
fileName = ['Solution map plot nRW = ' num2str(nRW) ' and np = ' num2str(fwd.np) ' and acc = ' num2str(acc) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);

%% Consider hte C_Bodin and D_Bodin
%% load the data
fwddataName = ['fwd_' num2str(k) 'k_type' num2str(fwd.mtype)];
load(fwddataName);
LABdataName = ['LAB_RW_Bodin_Result_' num2str(k) 'k_type' num2str(fwd.mtype)];
load(LABdataName);
%% Calcualte the stable chains
nburn = 0.2*nRW_Bodin;
t = 5;
C_stable_Bodin = C_Bodin(nburn:t:end,:);
D_stable_Bodin = D_Bodin(nburn:t:end,:);

%% Plot the LAB tomography for each sample
f = figure
for i = 1:5
    LABmesh = setParameterization(C_stable_Bodin(i,:),D_stable_Bodin(i,:),fwd);
    subplot(1,6,i)
    plotLAB(LABmesh,'plotSites',1);
    hold on
    plotMesh(fwd.mesh.X,fwd.mesh.T,'elementType',1);
    title('Sample ',i)
    xlabel('x');ylabel('z');
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis equal tight
end
subplot(1,6,6)
title('Average')
solutionMap(C_stable_Bodin,D_stable_Bodin,fwd)
x0=10;
y0=10;
width=1000;
height=350;
set(gcf,'position',[x0,y0,width,height])
fileName = ['LAB plot for first five samples with average solution Bodin nRW = ' num2str(nRW) ' and np = ' num2str(fwd.np) ' and acc = ' num2str(acc_Bodin) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
%% Plot the average solution map
solutionMap(C_stable_Bodin,D_stable_Bodin,fwd)
fileName = ['Solution map plot Bodin nRW = ' num2str(nRW_Bodin) ' and np = ' num2str(fwd.np) ' and acc = ' num2str(acc_Bodin) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
%% plot and average the image
% for i = 1:50
%     f = figure;
%     LABmesh = setParameterization(C_stable(i,:),D_stable(i,:),fwd);
%     plotLAB(LABmesh,'plotSites',1);
%     title('Sample ',i)
%     xlabel('x');ylabel('z');
%     set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
%     axis equal tight
%     fileName = ['LAB plot Bodin Sample = ' num2str(i) '.png'];
%     fn = fullfile(figurepath, fileName);
%     saveas(f,fn);
% end
% %%
% I0 = imread('figuresVoronoi/LAB plot Bodin Sample = 1.png');
% sumImage = double(I0);
% fileName = ['LAB plot Bodin Sample = ' num2str(i) '.png'];
% fn = fullfile(figurepath, fileName);
% for i = 2:50 % Read in remaining images
%     rgbImage = imread(fn);
%     %export_fig(rgbImage)
%     sumImage = sumImage + double(rgbImage);
% end
% meanImage = sumImage/10;
% meanImage = uint8(meanImage);
% figure;
% %plotMesh(fwd.mesh.X,fwd.mesh.T,'elementType',1);
% %hold on
% % plotMesh(fwd.mesh.XP,fwd.mesh.TP,'elementType',1,'labelNodes', 1 );
% imshow(meanImage);title('Average');
% 
% 
% %% plot and average the image
% for i = 1:50
%     f = figure;
%     LABmesh = setParameterization(C_stable(i,:),D_stable(i,:),fwd);
% %    subplot(1,5,i)
%     plotLAB(LABmesh,'plotSites',1);
%     hold on
%     plotMesh(fwd.mesh.X,fwd.mesh.T,'elementType',1);
%     title('Sample ',i)
%     xlabel('x');ylabel('z');
%     set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
%     axis equal tight
%     fileName = ['LAB and mesh plot Bodin Sample = ' num2str(i) '.png'];
%     fn = fullfile(figurepath, fileName);
%     saveas(f,fn);
% end
% %%
% I0 = imread('figuresVoronoi/LAB and mesh plot Bodin Sample = 1.png');
% sumImage = double(I0);
% fileName = ['LAB plot Bodin Sample = ' num2str(i) '.png'];
% fn = fullfile(figurepath, fileName);
% for i = 2:50 % Read in remaining images
%     rgbImage = imread(fn);
%     %export_fig(rgbImage)
%     sumImage = sumImage + double(rgbImage);
% end
% meanImage = sumImage/10;
% meanImage = uint8(meanImage);
% figure;
% %plotLAB(true_LABmesh,'plotSites',1);
% %hold on
% imshow(meanImage);title('Average');
