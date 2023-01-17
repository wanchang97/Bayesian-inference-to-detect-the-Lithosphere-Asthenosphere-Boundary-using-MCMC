function B0 = makeQoIMatrix(sensor,mesh)
%
nOfObservations = length(sensor);
nDof = mesh.nOfVelDof + mesh.nOfPreDof;

B0 = zeros(nOfObservations,nDof);

for ob = 1:nOfObservations
   node = sensor(ob).ixNode;

   switch (sensor(ob).component)
      case 1 % ux
         ix = 2*node-1;
      case 2 % uy
         ix = 2*node;
      case 3 % pres
         ix = mesh.nOfVelDof + node;
   end

   B0(ob,ix) = 1;
end
