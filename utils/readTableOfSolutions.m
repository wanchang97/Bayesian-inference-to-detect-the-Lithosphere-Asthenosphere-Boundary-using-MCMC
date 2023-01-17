function table = readTableOfSolutions(fwd)
%
fprintf('Reading forward solutions from file...\n')

np = fwd.param.nOfParam;
nx = fwd.mesh.nx;
nz = fwd.mesh.nz;
nq = fwd.mesh.refElem.nGaussPoints;

fileName = sprintf('data_%d_%d_%d_%d.mat', np,nx,nz,nq);
folderName = 'storedData';
fn = fullfile(folderName, fileName);

S = load(fn);

table.grid = S.grid;
table.solVector = S.solVector;