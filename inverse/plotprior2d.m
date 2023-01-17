function plotprior2d(fwd,mu,cov,type,grid,figurepath)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Evaluate the E and V using the Middle Points rule. however, with the limit
of our model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inputs:
    fwd:           prior pdf(m)
    Covmatrix_obs:   
    plotInfo: 
    figurepath:
Outputs: 
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%Covmatrix_obs = diag(sigma_ob.^2);

fprintf("Starting Plot Likelihood \n");
[X,Y] = meshgrid(grid,grid);
X_vec = X(:);
Y_vec = Y(:);
P = priorPdf([X_vec'; Y_vec']',mu',cov,type)
P = reshape(P, size(X));

figure
s = surf(X,Y,P);
%,'FaceAlpha',0.5
% s.EdgeColor = 'none';


colormap jet
lighting gouraud
shading interp


xlabel("m1");
ylabel("m2");
zlabel("Prior")
title("Prior pdf")

% cmin = fwd.param.cmin;
% cmax = fwd.param.cmax;
% dmin = fwd.param.dmin;
% dmax = fwd.param.dmax;
mmin = fwd.param.mmin;
mmax = fwd.param.mmax;
xlim([mmin,mmax]);
ylim([mmin,mmax]);
fileName = ["Prior_2d plot .png"]
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);