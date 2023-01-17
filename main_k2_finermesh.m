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
FE.nx = 100; FE.nz = 100;
% number of Gauss points in each element 3x3
FE.nGaussPoints = 9; 
fwd.mesh = createMeshForVelocityandPressure(fwd.geometry,FE);
% parametric domains       
fwd.param.dmin = 0.1 * fwd.geometry.lz;
fwd.param.dmax = 0.9 * fwd.geometry.lz;
fwd.param.bmin = 0;
fwd.param.bmax = fwd.geometry.lx;
%% Observables
sensorExactLocations = ...
    [fwd.geometry.lx/3,   fwd.geometry.lz*2/5  ; ...
    fwd.geometry.lx*1/4, fwd.geometry.lz*3/5  ; ...
    fwd.geometry.lx*1/5, fwd.geometry.lz*11/20; ...
    fwd.geometry.lx/4,   fwd.geometry.lz*9/20  ; ...
    fwd.geometry.lx/2,   fwd.geometry.lz*7/20  ];
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
k = 2; fwd.k = k;
% Parameterization type
% 0: Fixed discretization,fixed dimension m = d : np = k
% 1: Updated discretization, fixed dimension m = (B,d):the boundary between the columns: np = k +(nsd-1)*(k-1)
fwd.mtype = 1;
switch  fwd.mtype
    case 0
        fwd.np = k;
        true_d = fwd.geometry.lz*[1/3;1/2];
        true_m = true_d;
        figurepath = ['figures_k=' num2str(k) '_finermesh/regularGrids'];
        stringinput = {};
        for i = 1:k
            stringinput{end+1} = ['d' num2str(i)];
        end
        plotmin = ones(k,1)*fwd.param.dmin;
        plotmax = ones(k,1)*fwd.param.dmax;
    case 1
        fwd.np = k+(k-1)*(fwd.nsd-1);
        true_b = fwd.geometry.lx*3/8;
        true_d = fwd.geometry.lz*[1/3;1/2];
        true_m = [true_b; true_d];
        figurepath = ['figures_k=' num2str(k) '_finermesh/irregularGrids'];
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
plotLAB_modelerror(fwd,true_LABmesh,delta)
for I = 1:length(sensor)
   x = sensor(I).nodeLocation;
   if sensor(I).component == 1
       mark = '_r';
   else
       mark = '|r';
   end
   plot(x(1),x(2),mark,'markersize', 11, 'linewidth', 3);
end
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title("True LAB associated with FEM mesh")
fileName = ['True LAB set up, mesh plot and sensor location with np = ' num2str(fwd.np) ' and mtype =' num2str(fwd.mtype) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);
%% Inverse
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

%LABmesh_prior = setParameterization(mean_prior_m,fwd);
%sol_prior = forward(fwd,LABmesh_prior);
% test priorPdf
logpiPdf_test = logpiPdf(true_m);

% Likelihood
sigmaObratio = 0.2;
sigma_ob = abs(true_data) * sigmaObratio;
Covmatrix_obs = diag(sigma_ob.^2);
my_likelihood_fixsigma = @(m) likelihood(fwd,m,Covmatrix_obs);
logmy_likelihood_fixsigma = @(m) loglikelihood(fwd,m,Covmatrix_obs);
% test error and likelihood
Errtest_true = Error(fwd,true_m);

Ltest_true_1 = likelihood(fwd,true_m,Covmatrix_obs);
Ltest_true_2 = my_likelihood_fixsigma(true_m);

logLtest_truem_1 = loglikelihood(fwd,true_m,Covmatrix_obs);
logLtest_truem_2 = logmy_likelihood_fixsigma(true_m);
logLtest_prior = logmy_likelihood_fixsigma(mean_prior_m);

%% Random Walk sampler
nRW = 5000; % first 1000 burn in and second 1000
fprintf('\nStarting Random Walk Sampler approximation\n ');
rng('default')
sigmaRWratio = 1/200;
switch fwd.mtype
    case 0
        d0 = ones(k,1)*fwd.geometry.lz/2;
        m0 = d0;
        % cells coordinates and the depth should match each other /500 = 0.79
        sigmad = fwd.geometry.lz*sigmaRWratio ; % acc/nRW = 0.35 % plot the chains
        sigmam.d = sigmad;
    case 1
        b0 = ones(k-1,1)*fwd.geometry.lx/k;
        d0 = ones(k,1)*fwd.geometry.lz/2;
        m0 = [b0;d0];
        % cells coordinates and the depth should match each other /500 = 0.795
        sigmad = fwd.geometry.lz*sigmaRWratio ; % acc/nRW = 0.35 % plot the chains
        sigmab = fwd.geometry.lx*sigmaRWratio ;
        sigmam.b = sigmab;
        sigmam.d = sigmad;
end
% return all the generated samples of the Markov chain
DataName = ['LAB_RW_finermesh_Result_k=' num2str(k) '_mtype=' num2str(fwd.mtype) '_sigmaRratio=' num2str(sigmaRratio) '_sigmaObratio=' num2str(sigmaObratio) '_sigmaRWratio=' num2str(sigmaRWratio,'%.4f') '_nsamples=' num2str(nRW) '_log.mat'];

if isfile(DataName)
    load(DataName);
    fprintf(['\nLoad file ' DataName ' finished \n']);
else
    [MRW,lnLRW,accRW] = logRW(logpiPdf,logmy_likelihood_fixsigma,m0,nRW,sigmam,fwd);
    save(DataName,'MRW','lnLRW','nRW','accRW')
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
fileName_TraceHist = ['Histogram and trace plot of with np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio,'%.4f')  ' and n = ' num2str(nRW) '.png'];
TraceHistogramPlot(MRW,true_m,delta,m0,figurepath,fwd,fileName_TraceHist)
%% Trace plot
fileName_Trace = ['Markov Chain states plot np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio,'%.4f') ' and n = ' num2str(nRW) '.png'];
TracePlot(MRW,true_m,delta,m0,nRW,m_RW,accRW,figurepath,fwd,sigmam,fileName_Trace);

%% Histogram plot of all inputs
fileName_Hist = ['normalized histogram plot of RW samples without burn-in period np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio,'%.4f')  ' n = ' num2str(nRW) '.png'];
HistogramPlot(MRW,true_m,delta,error_m_RW,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,sigmam,sigmaObratio,fileName_Hist)

%% True LAB setup and the approximated LAB setup
figure,clf
subplot(1,2,1)
plotLAB_modelerror(fwd,true_LABmesh,delta)
% plotLAB(true_LABmesh)
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
title("true LAB")
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
title("recovered LAB using RW")
fileName = ['Comparison of LAB np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio)  ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio,'%.4f') ' and nRW = ' num2str(nRW) '.png'];
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);

%% Plot the histogram and the scatter plot together
fileName = ['Scatter plot of steady MCMC samples (without burn-in periods np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio,'%.4f') ' and nRW = ' num2str(nRW) '.png'];
ScatterPlot(MRW,true_m,nburn,nRW,m_RW,accRW,figurepath,fwd,stringinput,plotmin,plotmax,sigmam,sigmaObratio,fileName);
%% Plot the lnLikelihood traceplot and histogram plot
fileName = ['Trace Histogram plot of log likelihood all accepted samples np = ' num2str(fwd.np) ' sigmaRratio = ' num2str(sigmaRratio) ' sigmaObratio = ' num2str(sigmaObratio) '_sigmaRWratio = ' num2str(sigmaRWratio,'%.4f') ' and nRW = ' num2str(nRW) '.png'];
TraceHistogramlnLPlot(lnLRW(nburn:end),figurepath,fileName)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Adaptive Metropolis algorithm
% nAM = 2000;
% fprintf('\nStarting Adaptive Metropolis algorithm approximation\n ');
% rng('default')
% switch fwd.mtype
%     case 0
%         d0 = ones(k,1)*fwd.geometry.lz/2;
%         m0 = d0;
%         % cells coordinates and the depth should match each other /500 = 0.795
%         sigmad = fwd.geometry.lz/200 ; % acc/nRW = 0.35 % plot the chains
%         sigmam.d = sigmad;
%     case 1
%         b0 = ones(k-1,1)*fwd.geometry.lx/2;
%         d0 = ones(k,1)*fwd.geometry.lz/2;
%         m0 = [b0;d0];
%         % cells coordinates and the depth should match each other /500 = 0.795
%         sigmad = fwd.geometry.lz/200 ; % acc/nRW = 0.35 % plot the chains
%         sigmab = fwd.geometry.lx/200 ;
%         sigmam.b = sigmab;
%         sigmam.d = sigmad;
% end
% DataName = ['LAB_AM_Result_k=' num2str(k) '_mtype=' num2str(fwd.mtype) '_nsamples=' num2str(nAM) '_log.mat'];
% if isfile(DataName)
%     load(DataName);
% else
%     [MAM,accAM] = logAM(logpiPdf,logmy_likelihood_fixsigma,m0,nAM,sigmam,fwd);
%     save(DataName,'MAM','nAM','accAM')
% end
% BAM = MAM(:,1:k-1);
% DAM = MAM(:,k:end);
% %% evaluate the ensemble solution
% nburn = 0.5*nAM;
% [E_AM,Cov_AM,E_uncertainty_AM] = Ensemble_RW(MAM,fwd,nburn);
% b_AM = E_AM(1:k-1);
% d_AM = E_AM(k:end);
% LABmesh_AM = setParameterization(E_AM,fwd);
% obs_AM = fwd.B0*forward(fwd,LABmesh_AM);
% switch fwd.mtype
%     case 0
%         LABmesh_AM = setParameterization(E_AM,fwd);
%         obs_AM = fwd.B0*forward(fwd,LABmesh_AM);
%         fprintf('\nExpectation_AM: %g %g %g %g\n',E_AM)
%         fprintf('Relative errors_E: %.4f%% %.4f%% %.4f%% %.4f%%\n',abs((E_AM-true_m) ./ true_m)*100)
%     case 1
%         b_AM = E_AM(1:k-1);
%         d_AM = E_AM(k:end);
%         LABmesh_AM = setParameterization(E_AM,fwd);
%         obs_AM = fwd.B0*forward(fwd,LABmesh_AM);
%         fprintf('\nb_AM: %g %g %g %g\n',b_AM);
%         fprintf('\nd_AM: %g %g %g %g\n',d_AM);
%         fprintf('Relative errors_b: %.4f%% %.4f%% %.4f%% %.4f%%\n',abs((b_AM-true_b) ./ true_b)*100)
%         fprintf('Relative errors_d: %.4f%% %.4f%% %.4f%% %.4f%%\n',abs((d_AM-true_d) ./ true_d)*100)
% end
% 
% %% Trace plot and the histogram plot
% for i = 1:fwd.np
%     figure
%     subplot(1,2,1) % histogram plot
%     histogram(MAM(:,i),'Normalization','probability');
%     hold on
%     plot(true_m(i),0,'*','MarkerSize',20)
%     xlabel(['m_',num2str(i)])
%     ylabel(['\pi(m_',num2str(i),')'])
%     title('histogram plot')
%     camroll(90)
% 
%     subplot(1,2,2) % trace plot
%     plot(MAM(:,i));
%     hold on
%     yline(true_m(i),'r','LineStyle','-.');
%     xlabel('steps');
%     ylabel(['m',num2str(i)]);
%     title('trace plot')
%     fileName = ['AM histogram of m_', num2str(i),' plot and trace plot nAM = ' num2str(nAM) ' and np = ' num2str(fwd.np) '.png'];
%     fn = fullfile(figurepath, fileName);
%     saveas(gcf,fn);
% end
% %% trace plot of all inputs
% fig = figure;
% for i = 1:k-1
%     subplot(2,k,i)
%     plot(BAM(:,i));
%     hold on
%     plot(0,b0(i),'ro');
%     xlabel('steps')
%     ylabel(['b',num2str(i)]);
%     yline(true_b(i),'--');
% end
% for i = 1:k
%     subplot(2,k,k+i)
%     plot(DAM(:,i));
%     hold on
%     plot(0,d0(i),'ro');
%     ylabel(['d',num2str(i)]);
%     xlabel('steps')
%     yline(true_d(i),'--');
% end
% % add relative errors
% subplot(2,k,k)
% for i = 1: k-1
%     txtb =['Error_{b_' num2str(i) '} = ',num2str(abs((b_AM(i)-true_b(i)) ./ true_b(i)))];
%     text(1/(k-1)*(i-1),0.8,txtb)
% end
% for i = 1:k
%     txtd = ['Error_{d_' num2str(i) '} = ',num2str(abs((d_AM(i)-true_d(i)) ./ true_d(i)))];
%     text(1/(k)*(i-1),0.6,txtd)
% end
% axis('off')
% % Common xlabel
% han = axes(fig,'visible','off');
% han.Title.Visible = 'on';
% han.XLabel.Visible = 'on';
% %xlabel(han,'steps')
% title(han,'Adaptive Metropolis Algorith trace plot of all inputs')
% fileName = ['Adaptive Metropolis Algorith trace plot of all inputs with nAM = ' num2str(nAM) ' and np = ' num2str(fwd.np) '.png'];
% fn = fullfile(figurepath, fileName);
% saveas(gcf,fn);
% %% histogram plot of all inputs
% fig = figure;
% for i = 1:k-1
%     subplot(2,k,i)
%     histogram(BAM(:,i),'Normalization','probability');
%     hold on
%     title(['histogram of b',num2str(i)])
%     plot(true_b(i,1),0,'rx')
%     xlim([fwd.param.bmin,fwd.param.bmax])
% end
% 
% for i = 1:k
%     subplot(2,k,k+i)
%     histogram(DAM(:,i),'Normalization','probability');
%     hold on
%     title(['histogram of D',num2str(i)])
%     plot(true_d(i,1),0,'rx')
%     xlim([fwd.param.dmin,fwd.param.dmax])
% end
% han = axes(fig,'visible','off');
% han.Title.Visible = 'on';
% han.XLabel.Visible = 'on';
% title(han,'Normalized histogram plot of RW samples','Position',[0.5, -0.1, 0])
% fig.Position = [100 100 600 600];
% fileName = ['Normalized histogram plot of AM samples nAM = ' num2str(nAM) ' np = ' num2str(fwd.np) ' and mtype = ' num2str(fwd.mtype) '.png'];
% fn = fullfile(figurepath, fileName);
% saveas(gcf,fn);
% 
% 
% %% Comparison of LAB
% figure,clf
% subplot(1,2,1)
% plotLAB(true_LABmesh)
% hold on
% plotVelocityField(true_LABmesh,fwd)
% for I = 1:length(sensor)
%     x = sensor(I).nodeLocation;
%     if sensor(I).component == 1
%         mark = '_r';
%     else
%         mark = '|r';
%     end
%     plot(x(1),x(2),mark,'markersize',11,'linewidth',3);
% end
% title("true LAB")
% subplot(1,2,2)
% plotLABuncertainty(E_AM,E_uncertainty_AM,fwd)
% %plotLAB(LABmesh_AM);
% plotVelocityField(LABmesh_AM,fwd)
% for i = 1: k-1
%     txtb =['Error_{b_' num2str(i) '} = ',num2str(abs((b_AM(i)-true_b(i)) ./ true_b(i)))];
%     text(sum(b_AM(1:i))*0.9,(d_AM(i)+d_AM(i+1))/2,txtb)
% end
% for i = 1:k
%     txtd = ['Error_{d_' num2str(i) '} = ',num2str(abs((d_AM(i)-true_d(i)) ./ true_d(i)))];
%     text(sum(b_AM(1:i-1))*0.9,d_AM(i)*0.9,txtd)
% end
% title("recovered LAB using AM")
% fileName = ['Comparison of LAB nAM = ' num2str(nAM) ' and np = ' num2str(fwd.np) ' and mtype = ' num2str(fwd.mtype) '.png'];
% fn = fullfile(figurepath, fileName);
% saveas(gcf,fn);

