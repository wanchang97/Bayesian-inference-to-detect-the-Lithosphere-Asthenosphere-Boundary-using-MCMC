function table = createTableOfSolution(fwd,plotInfo)
%
fprintf('Saving forward solution to file...\n')

grid      = plotInfo.grid;
midPoints = (grid(1:end-1)+grid(2:end))/2;
nOfPoints = length(midPoints);

np = fwd.param.nOfParam;
nx = fwd.mesh.nx;
nz = fwd.mesh.nz;
nq = fwd.mesh.refElem.nGaussPoints;
%nob = size(fwd.observedData,1);

% solVector = zeros(nOfPoints^(np),nob);
% solVector (i,:) = forward(fwd,m())
% m = ()

a = nOfPoints*ones(1,np);
solVector = cell(a);
switch np
    case 3
        % FIXME: 3 parameters assumed how to write a general np case
        for j = 1 : nOfPoints
            for k = 1:nOfPoints
                for l = 1:nOfPoints
                    m = [midPoints(j) midPoints(k) midPoints(l)];
                    sol = forward(fwd, m);
                    solVector{j,k,l} = sol;
                end
            end
            fprintf('%d/%d\n',j,nOfPoints)
        end
    case 2
        % FIXME: 2 parameters assumed
        for p1 = 1:nOfPoints
            for p2 = 1:nOfPoints
                m = [midPoints(p1) midPoints(p2)];
                sol = forward(fwd, m);

                solVector{p1,p2} = sol;
            end
            fprintf('%d/%d\n',p1,nOfPoints)
        end
end



fileName = sprintf('data_%d_%d_%d_%d', np,nx,nz,nq);
folderName = 'storedData';
fn = fullfile(folderName, fileName);

save(fn, 'solVector', 'grid', 'midPoints','-v7.3')

table.grid = grid;
table.solVector = solVector;
