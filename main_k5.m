clearvars
close all
clc
%%
addpath("mesh");
addpath("utils");
addpath("inverse");
addpath("FE");
addpath("plot");
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
% parametric domains       
fwd.param.dmin = 0.1 * fwd.geometry.lz;
fwd.param.dmax = 0.9 * fwd.geometry.lz;
fwd.param.bmin = 0;
fwd.param.bmax = fwd.geometry.lx;
% plot the mesh 
figure,clf
plotMesh(fwd.mesh.X,fwd.mesh.T,'plotNodes',0,'lineSpecNodes','.k','labelNodes', 0 ,'labelNodesFontSize', 8,'labelNodesColor', 'k', ...
   'plotElements', 1,'lineSpecElements', '-c','lineThicknessElements', 0.5,'labelElements', 0,'labelElementsFontSize', 8,'labelElementsColor', 'b','elementType', 1);
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title("The Finite Element mesh set up")
fileName = ['The FE set up with nx = ' num2str(FE.nx)  ', nz = ' num2str(FE.nz) '.png'];
saveas(gcf,fileName);
%% Observables
sensorExactLocations = ...
    [fwd.geometry.lx*4/24,   fwd.geometry.lz*5/20  ; ...
    fwd.geometry.lx*8/24,    fwd.geometry.lz*5/20 ; ...
    fwd.geometry.lx*12/24,    fwd.geometry.lz*5/20; ...
    fwd.geometry.lx*16/24,   fwd.geometry.lz*5/20 ; ...
    fwd.geometry.lx*20/24,   fwd.geometry.lz*5/20; ...
    fwd.geometry.lx*4/24,   fwd.geometry.lz*7/20  ; ...
    fwd.geometry.lx*8/24,    fwd.geometry.lz*7/20 ; ...
    fwd.geometry.lx*12/24,    fwd.geometry.lz*7/20; ...
    fwd.geometry.lx*16/24,   fwd.geometry.lz*7/20 ; ...
    fwd.geometry.lx*20/24,   fwd.geometry.lz*7/20; ...
    fwd.geometry.lx*4/24,   fwd.geometry.lz*9/20  ; ...
    fwd.geometry.lx*8/24,    fwd.geometry.lz*9/20 ; ...
    fwd.geometry.lx*12/24,    fwd.geometry.lz*9/20; ...
    fwd.geometry.lx*16/24,   fwd.geometry.lz*9/20 ; ...
    fwd.geometry.lx*20/24,   fwd.geometry.lz*9/20; ...
    fwd.geometry.lx*4/24,   fwd.geometry.lz*11/20  ; ...
    fwd.geometry.lx*8/24,    fwd.geometry.lz*11/20 ; ...
    fwd.geometry.lx*12/24,    fwd.geometry.lz*11/20; ...
    fwd.geometry.lx*16/24,   fwd.geometry.lz*11/20 ; ...
    fwd.geometry.lx*20/24,   fwd.geometry.lz*11/20; ];

nSensorLocation = size(sensorExactLocations,1);

% Sensors (at the closest node)
for I = 1:size(sensorExactLocations,1)
    sensor(2*I-1) = setSensor(fwd.mesh,sensorExactLocations(I,:),1); % x
    sensor(2*I)   = setSensor(fwd.mesh,sensorExactLocations(I,:),2); % y
end

% Quantity of interest
% qoi_scaling_constant = 1E11;
fwd.B0 = makeQoIMatrix(sensor,fwd.mesh);

%% Given LAB 
% LAB discretization
k = 5; fwd.k = k;
% Parameterization type
% 0: Fixed discretization,fixed dimension m = d : nm = k
% 1: Updated discretization, fixed dimension m = (B,d):the boundary between the columns: nm = k +(nsd-1)*(k-1)
fwd.mtype = 1;
switch  fwd.mtype
    case 0
        true_d = fwd.geometry.lz*[1/3;1/2;2/5;1/4;1/3];   
        %dz1 = sqrt(0.6)/(2*FE.nz)*fwd.geometry.lz;% model error associated with Gauss points and the position of the true LAB
        %delta = [dz1,dz1];
        fwd.nm = k;
        true_m = true_d;
        figurepath = ['figures_k=' num2str(k) '/regularGrids'];
        stringinput = {};
        for i = 1:k
            stringinput{end+1} = ['d' num2str(i)];
        end
        plotmin = ones(k,1)*fwd.param.dmin;
        plotmax = ones(k,1)*fwd.param.dmax;

    case 1
        true_b = fwd.geometry.lx*[3/18;2/18;4/18;4/18];
        true_d = fwd.geometry.lz*[1/3;1/2;2/5;1/4;1/3];   
        fwd.nm = k+(k-1)*(fwd.nsd-1);
        true_m = [true_b; true_d];
        figurepath = ['figures_k=' num2str(k) '/irregularGrids'];
        stringinput = {};
        for i = 1:k-1
            stringinput{end+1} = ['b' num2str(i)];
        end
        for i = 1:k
            stringinput{end+1} = ['d' num2str(i)];
        end
        plotmin = [ones(k-1,1)*fwd.param.bmin;ones(k,1)*fwd.param.dmin];
        plotmax = [ones(k-1,1)*fwd.param.bmax;ones(k,1)*fwd.param.dmax];
end
delta = ModelResolution(FE,fwd,true_m);
true_LABmesh = setParameterization(true_m,fwd);
sol = forward(fwd,true_LABmesh);
true_data = fwd.B0*sol;
sigmaRratio = 0;
fwd.observedData = true_data; % no noise at the moment
%%
% plot the mesh and forward LAB setup 
figure,clf
plotMesh(fwd.mesh.X,fwd.mesh.T,'plotNodes',0,'lineSpecNodes','.k','labelNodes', 0 ,'labelNodesFontSize', 8,'labelNodesColor', 'k', ...
   'plotElements', 1,'lineSpecElements', '-c','lineThicknessElements', 0.5,'labelElements', 0,'labelElementsFontSize', 8,'labelElementsColor', 'b','elementType', 1);
hold on
plotLAB_modelerror(fwd,true_LABmesh,delta,0)
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title("The Finite Element mesh and the true LAB set up")
fileName = ['The FE and the true LAB set up with nx = ' num2str(FE.nx)  ', nz = ' num2str(FE.nz) '.png'];
saveas(gcf,fileName);
% plot the true set up to visulize the LAB and the sensor location
figure,clf
hold on
plotVelocityField(true_LABmesh,fwd)
plotLAB_modelerror(fwd,true_LABmesh,delta,0)
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title("True LAB setup and its velocity field")
fileName = ['True LAB set up and velocity with nm = ' num2str(fwd.nm) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
% plot the true set up to visulize the LAB and the sensor location
figure,clf
hold on
for I = 1:length(sensor)
    x = sensor(I).nodeLocation;
    if sensor(I).component == 1
        mark = '_r';
    else
        mark = '|r';
    end
    plot(x(1),x(2),mark,'markersize',11,'linewidth',3);
end
plotVelocityField(true_LABmesh,fwd)
plotLAB_modelerror(fwd,true_LABmesh,delta,0)
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title("True LAB associated with FEM mesh error")
fileName = ['True LAB set up velocity field and sensor location with nm = ' num2str(fwd.nm) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
%% Inverse problem
% Prior
switch fwd.mtype
    case 0
        priorInfo.dtype = 'uniform';
        d_dist = Dist(priorInfo.dtype ,'PAR',[fwd.param.dmin,fwd.param.dmax]);
        priorInfo.mdist = d_dist;
        mean_prior_m = d_dist.mean*ones(k,1);
    case 1
        priorInfo.btype = 'uniform';
        b_dist = Dist(priorInfo.btype ,'PAR',[fwd.param.bmin,fwd.param.bmax]);
        priorInfo.dtype = 'uniform';
        d_dist = Dist(priorInfo.dtype ,'PAR',[fwd.param.dmin,fwd.param.dmax]);
        priorInfo.mdist = [b_dist,d_dist];

        mean_prior_b = b_dist.mean*ones(k-1,1);
        mean_prior_d = d_dist.mean*ones(k,1);
        mean_prior_m = [mean_prior_b;mean_prior_d];
end
piPdf = @(m)priorPdf(m,priorInfo.mdist,fwd);
logpiPdf = @(m)logpriorPdf(m,priorInfo.mdist,fwd);

priorInfo.piPdf = piPdf;
priorInfo.logpiPdf = logpiPdf;

logpiPdf_test = logpiPdf(true_m);
%% Plot the prior marginal distribution

% Likelihood
sigmaObratio = 0.5;%1;
sigma_ob = abs(true_data) * sigmaObratio;
Covmatrix_obs = diag(sigma_ob.^2);
my_likelihood_fixsigma = @(m) likelihood(fwd,m,Covmatrix_obs);
logmy_likelihood_fixsigma = @(m) loglikelihood(fwd,m,Covmatrix_obs);
% test error and likelihood
Err_test = Error(fwd,true_m);

Ltest_true_1 = likelihood(fwd,true_m,Covmatrix_obs);
Ltest_true_2 = my_likelihood_fixsigma(true_m);

logLtest_truem_1 = loglikelihood(fwd,true_m,Covmatrix_obs);
logLtest_truem_2 = logmy_likelihood_fixsigma(true_m);
logLtest_prior = logmy_likelihood_fixsigma(mean_prior_m);

%% Random Walk sampler
nRW = 20000;
%nRW = 10000; % first 1000 burn in and second 1000
fprintf('\nStarting Random Walk Sampler written in logscale approximation\n ');
rng('default')
%sigmaRWratio = 1/100;1/150;1/120;1/110;1/105;1/200;1/220;1/240
%sigmaRWratio_b = 1/100; sigmaRWratio_d = 1/200;

%sigmaRWratio_b = 1/220; sigmaRWratio_d = 1/220;

sigmaRWratio_b = 1/150; sigmaRWratio_d = 1/90;
%sigmaRWratio_d = 1/100; 1/150; 1/120; 1/100;
switch fwd.mtype
    case 0
        d0 = ones(k,1)*fwd.geometry.lz/2;
        m0 = d0;
        % cells coordinates and the depth should match each other /500 = 0.795        
        sigmad = fwd.geometry.lz*sigmaRWratio_d ; % acc/nRW = 0.35 % plot the chains
        sigmam.d = sigmad;
        sigmaRWratiotext = ['_sigmaRWratio=' num2str(sigmaRWratio_d,'%.4f')];
        sigmaRWratioPlottext = [' sigmaRWratio =' num2str(sigmaRWratio_d,'%.4f')];
    case 1
        b0 = ones(k-1,1)*fwd.geometry.lx/k;
        d0 = ones(k,1)*fwd.geometry.lz/2;
        m0 = [b0;d0];
        logtest_0 = logmy_likelihood_fixsigma(m0);
        % cells coordinates and the depth should match each other /500 = 0.795
        sigmad = fwd.geometry.lz*sigmaRWratio_d ; % acc/nRW = 0.35 % plot the chains
        sigmab = fwd.geometry.lx*sigmaRWratio_b ;
        sigmam.b = sigmab;
        sigmam.d = sigmad;
        sigmaRWratiotext = ['_sigmaRWratio_b=' num2str(sigmaRWratio_b,'%.4f') '_sigmaRWratio_d=' num2str(sigmaRWratio_d,'%.4f')];
        sigmaRWratioPlottext = [' sigmaRWratio_b = ' num2str(sigmaRWratio_b,'%.4f') ' sigmaRWratio_d = ' num2str(sigmaRWratio_d,'%.4f')];
end
% return all the generated samples of the Markov chain
%DataName = ['LAB_RW_Result_k=' num2str(k) '_mtype=' num2str(fwd.mtype) '_sigmaRratio=' num2str(sigmaRratio) '_sigmaObratio=' num2str(sigmaObratio) '_sigmaRWratio_b=' num2str(sigmaRWratio_b,'%.4f') '_sigmaRWratio_d=' num2str(sigmaRWratio_d,'%.4f') '_nsamples=' num2str(nRW) '_log.mat'];
DataName = ['LAB_RW_Result_k=' num2str(k) '_mtype=' num2str(fwd.mtype) '_sigmaRratio=' num2str(sigmaRratio) '_sigmaObratio=' num2str(sigmaObratio)  sigmaRWratiotext  '_nsamples=' num2str(nRW) '_log.mat'];

if  isfile(DataName)
    load(DataName);
    fprintf(['\nLoad file ' DataName ' finished \n']);
else
    [MRW,lnLRW,accRW] = logRW(logpiPdf,logmy_likelihood_fixsigma,m0,nRW,sigmam,fwd);
    save(DataName,'MRW','lnLRW','nRW','accRW','sensor')
end
%% evaluate the ensemble solution
nburn = 0.5*nRW;
t = 1; % number of intermediate samples
alpha = 0.05 ;% 95% CI
[m_RW,Cov_RW,CI_RW] = Ensemble_RW(MRW,fwd,nburn,t,alpha);
% The autocorrelation matrix
temp = sqrt(diag(Cov_RW));
Corr_RW = Cov_RW./(temp*temp');
error_m_RW = abs((m_RW-true_m) ./ true_m)*100;
modelerror_m = (delta(2,:)+delta(1,:))'./m_RW;

switch fwd.mtype
    case 0
        plotmin2 = m_RW - 3*sqrt(diag(Cov_RW));
        plotmax2 = m_RW + 3*sqrt(diag(Cov_RW));
        fprintf('\nd_RW: %g %g %g %g\n',m_RW)
        fprintf('\nRelative errors_d: %.4f%% %.4f%% %.4f%% %.4f%%\n',error_m_RW)
        fprintf('\nRelative FEM errors_d: %.4f%% %.4f%% %.4f%% %.4f%% \n',modelerror_m)
    case 1
        plotmin2 = m_RW - 3*sqrt(diag(Cov_RW));
        plotmax2 = m_RW + 3*sqrt(diag(Cov_RW));
        BRW = MRW(:,1:k-1);
        DRW = MRW(:,k:end);
        b_RW = m_RW(1:k-1);
        d_RW = m_RW(k:end);
        error_b_RW = error_m_RW(1:k-1);
        error_d_RW = error_m_RW(k:end);

        fprintf('\nb_RW: %g %g %g %g\n',b_RW);
        fprintf('\nd_RW: %g %g %g %g\n',d_RW);
        fprintf('\nRelative errors_b: %.4f%% %.4f%% %.4f%% %.4f%%\n',error_b_RW)
        fprintf('\nRelative errors_d: %.4f%% %.4f%% %.4f%% %.4f%%\n',error_d_RW)
        fprintf('\nRelative FEM errors_b: %.4f%% %.4f%% %.4f%% %.4f%% \n',modelerror_m(1:k-1));
        fprintf('\nRelative FEM errors_d: %.4f%% %.4f%% %.4f%% %.4f%% \n',modelerror_m(k:end));
end
% Predictive solution using point estimator
LABmesh_RW = setParameterization(m_RW,fwd);
obs_RW = fwd.B0*forward(fwd,LABmesh_RW);
fprintf('\nNorm of the difference between the predictive and the true observations  %.4f%%  \n',norm(obs_RW-true_data)/norm(true_data)*100);

%% Trace plot and the histogram plot for all samples
%fileName_TraceHist = ['Histogram and trace plot with nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) ' sigmaRWratio_b = ' num2str(sigmaRWratio_b,'%.4f') ' sigmaRWratio_d = ' num2str(sigmaRWratio_d,'%.4f')  ' and n = ' num2str(nRW) '.png'];
fileName_TraceHist = ['Histogram and trace plot with nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext ' and n = ' num2str(nRW) '.png'];

TraceHistogramPlot(MRW,true_m,delta,m0,figurepath,fwd,fileName_TraceHist)
%% Trace plot for all samples
fileName_Trace = ['Markov Chain states plot nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObRatio = ' num2str(sigmaObratio) sigmaRWratioPlottext ' and n = ' num2str(nRW) '.png'];
TracePlot(MRW,true_m,delta,m0,nRW,m_RW,accRW,figurepath,fwd,sigmam,fileName_Trace);

%% Histogram plot of all inputs without burn-in periods
fileName_Hist = ['normalized histogram plot of RW samples without burn-in period with nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext  'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn) ' n = ' num2str(nRW) '.png'];
HistogramPlot(MRW,true_m,delta,error_m_RW,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName_Hist,priorInfo)

fileName_Hist = ['normalized histogram plot of RW samples without burn-in period with prior pdf nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext  'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn) ' n = ' num2str(nRW) '.png'];
HistogramPlot0(MRW,true_m,delta,error_m_RW,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName_Hist,priorInfo)

fileName_Hist = ['normalized histfit plot of RW samples without burn-in period with prior pdf nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext 'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn) ' n = ' num2str(nRW) '.png'];
HistogramPlot1(MRW,true_m,delta,error_m_RW,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName_Hist,priorInfo)

fileName_Hist = ['normalized histfit plot of RW samples for burn-in period nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext 'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn) ' n = ' num2str(nRW) '.png'];
HistogramPlot2(MRW,true_m,delta,error_m_RW,nburn/10,nRW,m_RW,accRW,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName_Hist)
%% True LAB setup and the approximated LAB setup
figure,clf
subplot(1,2,1)
plotLAB_modelerror(fwd,true_LABmesh,delta,0)
hold on
plotVelocityField(true_LABmesh,fwd)
for I = 1:length(sensor)
    x = sensor(I).nodeLocation;
    if sensor(I).component == 1
        mark = '_r';
    else
        mark = '|r';
    end
    plot(x(1),x(2),mark,'markersize',11,'linewidth',3);
end
title("true LAB with FEM model error")
subplot(1,2,2)
plotLABuncertainty(m_RW,CI_RW,fwd)
plotVelocityField(LABmesh_RW,fwd)
title("estimated LAB and the Confidence Interval")
fileName = ['Comparison of LAB nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio)  ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext 'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn) ' and nRW = ' num2str(nRW) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);

 %%% Plot the histogram and the scatter plot together without burnin periods
 %fileName = ['Scatter plot of steady MCMC samples (without burn-in periods nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) ' sigmaRWratio_b = ' num2str(sigmaRWratio_b,'%.4f') ' sigmaRWratio_d = ' num2str(sigmaRWratio_d,'%.4f') 'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn) ' and nRW = ' num2str(nRW) '.png'];
 %ScatterPlot(MRW,true_m,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,plotmin,plotmax,sigmam,sigmaObratio,fileName);
%% Plot the histogram and the scatter plot together without burnin periods
fileName = ['Scatter plot of steady MCMC samples (without burn-in periods nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext 'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn) ' and nRW = ' num2str(nRW) '.png'];
t2 = 1;
ScatterPlot(MRW,nburn,t2,m_RW,figurepath,fwd,stringinput,plotmin2,plotmax2,fileName)
%% Plot the lnLikelihood traceplot and histogram plot without burnin periods
fileName = ['Trace Histogram plot of log likelihood all accepted samples nm = ' num2str(fwd.nm) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) sigmaRWratioPlottext 'every t = ' num2str(t) 'samples and nburn = ' num2str(nburn)  ' and nRW = ' num2str(nRW) '.png'];
TraceHistogramlnLPlot(lnLRW(nburn:end),figurepath,fileName)