function solutionMap(C_stable,D_stable,fwd)
% Asign all the depth D to our lx discretization according to C
% Average all depth on each pixel ( or FEM nodes)
nsamples = size(D_stable,1);
k = size(D_stable,2);
% depth of all pixels
pixelDepth = zeros(nsamples,fwd.mesh.nx);
% uniform x-coordinates
X_mid = (fwd.mesh.XP(1:fwd.mesh.nx,:)+fwd.mesh.XP(2:fwd.mesh.nx+1,:))/2;
locationx = X_mid(:,1)';
for i = 1:nsamples
    LABmesh = setParameterization(C_stable(i,:),D_stable(i,:),fwd);
    pixeld = LABdepth2(D_stable(i,:)',locationx,LABmesh);
    pixelDepth(i,:) = pixeld;
end
averagePixelDepth = mean(pixelDepth,1);
% plot the LAB depth based on the meshes

plotLABpixel(fwd.mesh.XP(1:fwd.mesh.nx+1,:),averagePixelDepth);
hold on
%plotLAB(fwd.true_LABmesh,'plotSites',0)
plotMesh(fwd.mesh.X,fwd.mesh.T,'elementType',1);
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
axis equal tight
end