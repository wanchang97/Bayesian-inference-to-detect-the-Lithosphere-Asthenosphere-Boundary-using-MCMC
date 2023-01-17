function [g,rho,nu] = getProperties(material,labDepth,z_gp)
%
idMaterial = (z_gp > labDepth') + 1;
g = material.gravity;
rho = material.rho(idMaterial);
nu = material.nu(idMaterial);
end

