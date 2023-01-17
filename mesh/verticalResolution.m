function gridPoints = verticalResolution(mesh)
%
% a square mesh is assumed
%
xmin = min(mesh.X(:,1));
xmax = max(mesh.X(:,1));
tol  = (xmax-xmin) / 1000;

% id of the nodes of the first column
nodes  = find(abs(mesh.XP(:,1)-xmin) < tol);
[~,ix] = sort(mesh.XP(nodes,2));
nodes  = nodes(ix);

n = round(sqrt(mesh.refElem.nGaussPoints));

zq = quadrature(0,n);
N  = shapeFunctions(0,2,zq);

gridPoints = mesh.XP(nodes(1),2);
for I = 1:length(nodes)-1
   Xe = [mesh.XP(nodes(I+1),2); mesh.XP(nodes(I),2)];
   gridPoints = [gridPoints; N*Xe];
end
gridPoints = [gridPoints; mesh.XP(nodes(end),2)];

% | m g  m g m  g  |    g    g    g  m1 | m2  g    g    g      |
