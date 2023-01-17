function plotLikelihood2d(fwd,Covmatrix_obs,grid,figurepath)
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
L = likelihood(fwd,[X_vec'; Y_vec'],Covmatrix_obs);
L = reshape(L, size(X));

figure
s = surf(X,Y,L);
%,'FaceAlpha',0.5
% s.EdgeColor = 'none';


colormap jet
lighting gouraud
shading interp


xlabel("m1");
ylabel("m2");
zlabel("Likelihood")
title("Likelihood pdf")

mmin = fwd.param.mmin;
mmax = fwd.param.mmax;

xlim([mmin,mmax]);
ylim([mmin,mmax]);
fileName = ["Likelihood_2d plot .png"]
%fileName = ["Likelihood_2d plot with sigma_ob = " num2str(sigma_ob) " .png"]
fn = fullfile(figurepath, fileName);
saveas(gcf,fn);