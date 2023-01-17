function LABmesh = setParameterization(m,fwd)
k = fwd.k;
xmin = 0;
xmax = fwd.geometry.lx;
switch fwd.mtype
    case 0 % 0: Fixed discretization,fixed dimension m = d : np = k
        LABmesh.T = [];
        LABmesh.b = (xmax-xmin)/k;% width
        LABmesh.X = linspace(xmin,xmax,k+1)';
        LABmesh.d = m;
    case 1 % Updated discretization, fixed dimension m = (b,d):the width and depth of columns: np = k +(nsd-1)*(k-1)
        b = m(1:k-1);
        d = m(k:end);
        BoundaryLoc = zeros(k-1,1);
        % for 1D
        for i = 1:k-1
            BoundaryLoc(i) = sum(b(1:i));
        end
        LABmesh.T = [];
        LABmesh.X = [xmin;BoundaryLoc;xmax];
        LABmesh.d = d;
        LABmesh.b = b;
end