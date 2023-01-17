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
%% Observables
% sensorExactLocations = ...
%     [fwd.geometry.lx*2/7,  fwd.geometry.lz*4/24  ; ...
%      fwd.geometry.lx*2/7,  fwd.geometry.lz*5/24  ; ...
%      fwd.geometry.lx*2/7,  fwd.geometry.lz*6/24  ; ...
%      fwd.geometry.lx*2/7,  fwd.geometry.lz*7/24  ; ...
%      fwd.geometry.lx*2/7,  fwd.geometry.lz*8/24  ; ...
%      fwd.geometry.lx*3/7,  fwd.geometry.lz*4/24  ; ...
%      fwd.geometry.lx*3/7,  fwd.geometry.lz*5/24  ; ...
%      fwd.geometry.lx*3/7,  fwd.geometry.lz*6/24  ; ...
%      fwd.geometry.lx*3/7,  fwd.geometry.lz*7/24  ; ...
%      fwd.geometry.lx*3/7,  fwd.geometry.lz*8/24  ; ...
%      fwd.geometry.lx*4/7,  fwd.geometry.lz*4/24  ; ...
%      fwd.geometry.lx*4/7,  fwd.geometry.lz*5/24  ; ...
%      fwd.geometry.lx*4/7,  fwd.geometry.lz*6/24  ; ...
%      fwd.geometry.lx*4/7,  fwd.geometry.lz*7/24  ; ...
%      fwd.geometry.lx*4/7,  fwd.geometry.lz*8/24  ; ...
%      fwd.geometry.lx*5/7,  fwd.geometry.lz*4/24  ; ...
%      fwd.geometry.lx*5/7,  fwd.geometry.lz*5/24  ; ...
%      fwd.geometry.lx*5/7,  fwd.geometry.lz*6/24  ; ...
%      fwd.geometry.lx*5/7,  fwd.geometry.lz*7/24  ; ...
%      fwd.geometry.lx*5/7,  fwd.geometry.lz*8/24  ; ...
%      fwd.geometry.lx*6/7,  fwd.geometry.lz*4/24  ; ...
%      fwd.geometry.lx*6/7,  fwd.geometry.lz*5/24  ; ...
%      fwd.geometry.lx*6/7,  fwd.geometry.lz*6/24  ; ...
%      fwd.geometry.lx*6/7,  fwd.geometry.lz*7/24  ; ...
%      fwd.geometry.lx*6/7,  fwd.geometry.lz*8/24];
% sensorExactLocations = ...
%     [fwd.geometry.lx*2/7,  fwd.geometry.lz*3/24  ; ...
%      fwd.geometry.lx*3/7,  fwd.geometry.lz*5/24  ; ...
%      fwd.geometry.lx*4/7,  fwd.geometry.lz*6/24  ; ...
%      fwd.geometry.lx*5/7,  fwd.geometry.lz*7/24  ; ...
%      fwd.geometry.lx*6/7,  fwd.geometry.lz*8/24];
% sensorExactLocations = ...
%     [fwd.geometry.lx*4/24,   fwd.geometry.lz*5/20 ; ...
%      fwd.geometry.lx*6/24,   fwd.geometry.lz*6/20 ; ...
%      fwd.geometry.lx*8/24,   fwd.geometry.lz*7/20 ; ...
%      fwd.geometry.lx*10/24,  fwd.geometry.lz*8/20 ; ...
%      fwd.geometry.lx*12/24,  fwd.geometry.lz*9/20 ; ...
%      fwd.geometry.lx*14/24,  fwd.geometry.lz*10/20; ...
%      fwd.geometry.lx*16/24,  fwd.geometry.lz*11/20; ...
%      fwd.geometry.lx*18/24,  fwd.geometry.lz*12/20; ...
%      fwd.geometry.lx*20/24,  fwd.geometry.lz*13/20];
% sensorExactLocations = ...
%     [fwd.geometry.lx*4/24,   fwd.geometry.lz*5/20  ; ...
%     fwd.geometry.lx*5/24,   fwd.geometry.lz*5.5/20  ; ...
%     fwd.geometry.lx*6/24,   fwd.geometry.lz*6/20  ; ...
%     fwd.geometry.lx*7/24,   fwd.geometry.lz*6.5/20  ; ...
%     fwd.geometry.lx*8/24,    fwd.geometry.lz*7/20 ; ...
%     fwd.geometry.lx*9/24,   fwd.geometry.lz*7.5/20  ; ...
%     fwd.geometry.lx*10/24,   fwd.geometry.lz*8/20  ; ...
%     fwd.geometry.lx*11/24,   fwd.geometry.lz*8.5/20  ; ...
%     fwd.geometry.lx*12/24,    fwd.geometry.lz*9/20; ...
%     fwd.geometry.lx*13/24,   fwd.geometry.lz*9.5/20  ; ...
%     fwd.geometry.lx*14/24,   fwd.geometry.lz*10/20  ; ...
%     fwd.geometry.lx*15/24,   fwd.geometry.lz*10.5/20  ; ...
%     fwd.geometry.lx*16/24,   fwd.geometry.lz*11/20 ; ...
%     fwd.geometry.lx*17/24,   fwd.geometry.lz*11.5/20  ; ...
%     fwd.geometry.lx*18/24,   fwd.geometry.lz*9/20  ; ...
%     fwd.geometry.lx*19/24,   fwd.geometry.lz*9.5/20  ; ...
%     fwd.geometry.lx*20/24,   fwd.geometry.lz*11/20; ];
% sensorExactLocations = ...
%     [fwd.geometry.lx*4/24,   fwd.geometry.lz*5/20  ; ...
%     fwd.geometry.lx*8/24,    fwd.geometry.lz*7/20 ; ...
%     fwd.geometry.lx*12/24,    fwd.geometry.lz*9/20; ...
%     fwd.geometry.lx*16/24,   fwd.geometry.lz*11/20 ; ...
%     fwd.geometry.lx*20/24,   fwd.geometry.lz*11/20; ];
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
% sensorExactLocations = ...
%     [fwd.geometry.lx*4/24,   fwd.geometry.lz*7/20  ; ...
%     fwd.geometry.lx*6/24,   fwd.geometry.lz*7/20  ; ...
%     fwd.geometry.lx*8/24,    fwd.geometry.lz*7/20 ; ...
%     fwd.geometry.lx*10/24,    fwd.geometry.lz*7/20 ; ...
%     fwd.geometry.lx*12/24,    fwd.geometry.lz*7/20; ...
%     fwd.geometry.lx*14/24,   fwd.geometry.lz*7/20; ...
%     fwd.geometry.lx*16/24,   fwd.geometry.lz*7/20 ; ...
%     fwd.geometry.lx*18/24,   fwd.geometry.lz*7/20; ...
%     fwd.geometry.lx*20/24,   fwd.geometry.lz*7/20; ...
%     fwd.geometry.lx*4/24,   fwd.geometry.lz*9/20  ; ...
%     fwd.geometry.lx*6/24,   fwd.geometry.lz*9/20  ; ...
%     fwd.geometry.lx*8/24,    fwd.geometry.lz*9/20 ; ...
%     fwd.geometry.lx*10/24,    fwd.geometry.lz*9/20 ; ...
%     fwd.geometry.lx*12/24,    fwd.geometry.lz*9/20; ...
%     fwd.geometry.lx*14/24,   fwd.geometry.lz*9/20; ...
%     fwd.geometry.lx*16/24,   fwd.geometry.lz*9/20 ; ...
%     fwd.geometry.lx*18/24,   fwd.geometry.lz*9/20; ...
%     fwd.geometry.lx*20/24,   fwd.geometry.lz*9/20; ...
%     fwd.geometry.lx*4/24,   fwd.geometry.lz*11/20  ; ...
%     fwd.geometry.lx*6/24,   fwd.geometry.lz*11/20  ; ...
%     fwd.geometry.lx*8/24,    fwd.geometry.lz*11/20 ; ...
%     fwd.geometry.lx*10/24,    fwd.geometry.lz*11/20 ; ...
%     fwd.geometry.lx*12/24,    fwd.geometry.lz*11/20; ...
%     fwd.geometry.lx*14/24,   fwd.geometry.lz*11/20; ...
%     fwd.geometry.lx*16/24,   fwd.geometry.lz*11/20 ; ...
%     fwd.geometry.lx*18/24,   fwd.geometry.lz*11/20; ...
%     fwd.geometry.lx*20/24,   fwd.geometry.lz*11/20; ];
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
k = 3; fwd.k = k;
% Parameterization type
% 0: Fixed discretization,fixed dimension m = d : np = k
% 1: Updated discretization, fixed dimension m = (B,d):the boundary between the columns: np = k +(nsd-1)*(k-1)
fwd.mtype = 0;
switch  fwd.mtype
    case 0
        true_d = fwd.geometry.lz*[1/3;1/2;2/5];   
        %dz1 = sqrt(0.6)/(2*FE.nz)*fwd.geometry.lz;% model error associated with Gauss points and the position of the true LAB
        %delta = [dz1,dz1];
        fwd.np = k;
        true_m = true_d;
        figurepath = ['figures_k=' num2str(k) '/regularGrids'];
        stringinput = {};
        for i = 1:k
            stringinput{end+1} = ['d' num2str(i)];
        end
        plotmin = ones(k,1)*fwd.param.dmin;
        plotmax = ones(k,1)*fwd.param.dmax;
    case 1
        true_b = fwd.geometry.lx*[3/8;3/12];
        true_d = fwd.geometry.lz*[1/3;1/2;2/5];   
        fwd.np = k+(k-1)*(fwd.nsd-1);
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

% plot the true set up to visulize the LAB and the sensor location
figure,clf
plotMesh(fwd.mesh.X,fwd.mesh.T,'plotNodes',0,'lineSpecNodes','.k','labelNodes', 0 ,'labelNodesFontSize', 8,'labelNodesColor', 'k', ...
   'plotElements', 1,'lineSpecElements', '-b','lineThicknessElements', 0.5,'labelElements', 0,'labelElementsFontSize', 8,'labelElementsColor', 'b','elementType', 1);
hold on
plotVelocityField(true_LABmesh,fwd)
plotLAB_modelerror(fwd,true_LABmesh,delta)
for I = 1:length(sensor)
    x = sensor(I).nodeLocation;
    if sensor(I).component == 1
        mark = '_r';
    else
        mark = '|r';
    end
    plot(x(1),x(2),mark,'markersize',11,'linewidth',3);
end
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title("True LAB associated with FEM mesh error")
fileName = ['True LAB set up, mesh plot and sensor location with np = ' num2str(fwd.np) '.png'];
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

% LABmesh_prior = setParameterization(mean_prior_m,fwd);
% sol_prior = forward(fwd,LABmesh_prior);
% test priorPdf
logpiPdf_test = logpiPdf(true_m);

% Likelihood
sigmaObratio = 0.5;
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
nRW = 10000; % first 1000 burn in and second 1000
fprintf('\nStarting Random Walk Sampler written in logscale approximation\n ');
rng('default')
%sigmaRWratio = 1/200;
sigmaRWratio = 1/120;
switch fwd.mtype
    case 0
        d0 = ones(k,1)*fwd.geometry.lz/2;
        m0 = d0;
        % cells coordinates and the depth should match each other /500 = 0.795        
        sigmad = fwd.geometry.lz*sigmaRWratio ; % acc/nRW = 0.35 % plot the chains
        sigmam.d = sigmad;
    case 1
        b0 = ones(k-1,1)*fwd.geometry.lx/k;
        d0 = ones(k,1)*fwd.geometry.lz/2;
        m0 = [b0;d0];
        logtest_0 = logmy_likelihood_fixsigma(m0);
        % cells coordinates and the depth should match each other /500 = 0.795
        sigmad = fwd.geometry.lz*sigmaRWratio ; % acc/nRW = 0.35 % plot the chains
        sigmab = fwd.geometry.lx*sigmaRWratio ;
        sigmam.b = sigmab;
        sigmam.d = sigmad;
end
% return all the generated samples of the Markov chain
% DataName = ['LAB_RW_Result_k=' num2str(k) '_mtype=' num2str(fwd.mtype) '_sigmaRratio=' num2str(sigmaRratio) '_sigmaObratio=' num2str(sigmaObratio) '_sigmaRWratio=' num2str(sigmaRWratio,'%.4f') '_nsamples=' num2str(nRW) '_log.mat'];
DataName = ['LAB_RW_Result_k=' num2str(k) '_mtype=' num2str(fwd.mtype) '_sigmaRratio=' num2str(sigmaRratio) '_sigmaObratio=' num2str(sigmaObratio) '_sigmaRWratio=' num2str(sigmaRWratio,'%.4f') '_nsamples=' num2str(nRW) '_log.mat'];

if isfile(DataName)
    load(DataName);
    fprintf(['\nLoad file ' DataName ' finished \n']);
else
    %[MallRW,MRW,accRW] = logRW_allsamples(logpiPdf,logmy_likelihood_fixsigma,m0,nRW,sigmam,fwd);
    [MRW,lnLRW,accRW] = logRW(logpiPdf,logmy_likelihood_fixsigma,m0,nRW,sigmam,fwd);
    %save(DataName,'MallRW','MRW','nRW','accRW')
    save(DataName,'MRW','lnLRW','nRW','accRW','sensor')
end

%% evaluate the ensemble solution
nburn = 0.5*nRW;
[m_RW,Cov_RW,E_uncertainty_RW] = Ensemble_RW(MRW,fwd,nburn);
error_m_RW = abs((m_RW-true_m) ./ true_m)*100;
switch fwd.mtype
    case 0
        fprintf('\nd_RW: %g %g %g %g\n',m_RW)
        fprintf('\nRelative errors_d: %.4f%% %.4f%% %.4f%% %.4f%%\n',error_m_RW)
    case 1
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
end
% Predictive solution using point estimator
LABmesh_RW = setParameterization(m_RW,fwd);
obs_RW = fwd.B0*forward(fwd,LABmesh_RW);
fprintf('\nNorm of the difference between the Predictive and the true observations  %.4f%%  \n',norm(obs_RW-true_data)/norm(true_data)*100);

%% Trace plot and the histogram plot
fileName_TraceHist = ['Histogram and trace plot with np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio)  ' and n = ' num2str(nRW) '.png'];
TraceHistogramPlot(MRW,true_m,delta,m0,figurepath,fwd,fileName_TraceHist)
%% Trace plot
fileName_Trace = ['Markov Chain states plot np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObRatio = ' num2str(sigmaObratio) ' sigmaRWatio = ' num2str(sigmaRWratio) ' and n = ' num2str(nRW) '.png'];
TracePlot(MRW,true_m,delta,m0,nRW,m_RW,accRW,figurepath,fwd,sigmam,fileName_Trace);

%% Histogram plot of all inputs
fileName_Hist = ['normalized histogram plot of RW samples without burn-in period np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio)  ' n = ' num2str(nRW) '.png'];
HistogramPlot(MRW,true_m,delta,error_m_RW,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName_Hist)

%% True LAB setup and the approximated LAB setup
figure,clf
subplot(1,2,1)
plotLAB_modelerror(fwd,true_LABmesh,delta)
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
plotLABuncertainty(m_RW,E_uncertainty_RW,fwd)
plotVelocityField(LABmesh_RW,fwd)
switch fwd.mtype
    case 0
        for i = 1:k
            txtd = ['Error_{d_' num2str(i) '} = ',num2str(error_m_RW(i),'%.4f%%')];
            text(fwd.geometry.lx/k*(i-1),m_RW(i)*0.9,txtd)
        end
    case 1
        for i = 1: k-1
            txtb =['Error_{b_' num2str(i) '} = ',num2str(error_b_RW(i),'%.4f%%')];
            text(sum(b_RW(1:i-1))*0.9,(d_RW(i)+fwd.geometry.lz)/2,txtb)
        end
        for i = 1:k
            txtd = ['Error_{d_' num2str(i) '} = ',num2str(error_d_RW(i),'%.4f%%')];
            text(sum(b_RW(1:i-1))*0.9,d_RW(i)*0.9,txtd)
        end
end
title("estimated LAB and the confidence interval using RW")
fileName = ['Comparison of LAB np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio)  ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio) ' and nRW = ' num2str(nRW) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);

%% Plot the histogram and the scatter plot together
fileName = ['Scatter plot of steady MCMC samples (without burn-in periods np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) ' sigmaRWratio = ' num2str(sigmaRWratio) ' and nRW = ' num2str(nRW) '.png'];
ScatterPlot(MRW,true_m,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,plotmin,plotmax,sigmam,sigmaObratio,fileName);

%% Plot the lnLikelihood traceplot and histogram plot
fileName = ['Trace Histogram plot of log likelihood all accepted samples np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) ' sigmaRWratio = ' num2str(sigmaRWratio,'%.4f') ' and nRW = ' num2str(nRW) '.png'];
TraceHistogramlnLPlot(lnLRW(nburn:end),figurepath,fileName)