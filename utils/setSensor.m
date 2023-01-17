function S = setSensor(mesh,location,component)
%

% Find closest node
n = size(mesh.X,1);
aux = mesh.X - repmat(location,n,1);
distanceToNodes = sqrt(sum(aux.^2,2));
[~,ix] = min(distanceToNodes);

S.ixNode = ix;
S.exactLocation = location;
S.nodeLocation = mesh.X(ix,:);
S.component = component;



