function plotLikelihood3d(fwd,covMatrix,plotInfo)

fprintf("Starting Plot Likelihood \n");
grid = plotInfo.grid;
nOfgrid = size(grid,1);
np = fwd.param.nOfParam;
m = zeros(nOfgrid^np,np);
for i = 1:nOfgrid
    for j = 1:nOfgrid
        for k = 1:nOfgrid
            m(nOfgrid^2*(i-1)+nOfgrid*(j-1)+k,:) = [grid(i) grid(j) grid(k)];
        end
    end
end
% m number of nsamples x np = nOfgrid^np x np
Li = likelihood(fwd,m,covMatrix);
L = reshape(Li,nOfgrid,nOfgrid,nOfgrid);
%figure
%plot(nOfgrid,L(:,))



