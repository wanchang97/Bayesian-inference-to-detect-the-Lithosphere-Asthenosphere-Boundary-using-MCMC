function plotLAB(labMesh)
X = labMesh.X;
d = labMesh.d;
n = length(X)-1;
%plotSites = getValueFromVarargin('plotSites',varargin,0);
hold on
for I = 1:n
    xx = [X(I) X(I+1)];
    yy = d(I) * [1 1];
    plot(xx,yy,'b','LineWidth',3);
end
end