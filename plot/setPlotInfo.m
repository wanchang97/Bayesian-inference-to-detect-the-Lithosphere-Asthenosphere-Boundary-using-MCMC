function plotInfo = setPlotInfo(fwd)

% grid of the parametric space
mm = verticalResolution(fwd.mesh);
plotInfo.grid = mm;
%
% FIXME: this only works for 2 paramenters
% [plotInfo.X,plotInfo.Y] = meshgrid(mm,mm);
% plotInfo.X_vec = plotInfo.X(:);
% plotInfo.Y_vec = plotInfo.Y(:);


