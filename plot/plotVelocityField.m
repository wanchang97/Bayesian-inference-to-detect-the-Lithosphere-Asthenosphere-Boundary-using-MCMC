function plotVelocityField(labMesh,fwd)
sol = forward(fwd,labMesh);
v = sol(1:fwd.mesh.nOfVelDof);
v = reshape(v',2,[])';
% Plot the velocity
quiver(fwd.mesh.X(:,1),fwd.mesh.X(:,2),v(:,1),v(:,2),"Color","#0072BD");
xlabel('x');ylabel('z');
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
axis equal tight
end
