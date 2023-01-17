function depth = LABdepth(locationx,labMesh)
%{
To find the LABdepth of all gauss points in the Finite element
Here the code is for the (Nearest Neighbor)Voronoin interpolation,discrete
For future we can have (linear)Delaunay interpolation or natural neighbor
(Voronois +  weighted average)
Input: 
    locationx ------ the given location | dim : (nsd-1) x nOfGauss 
    labMesh   .C    ------ LAB depth discretization |dim : k x (nsd-1)
              .d    ------ LAB depth discretization |dim : k x 1
              .T    ------ LAB depth discretization |dim : k x ()
              .X    ------ LAB depth discretization |dim : (k+1) x (nsd-1)

Output: 
    depth     ------ the LAB depth of all the Gausspoints in the element
                                         |dim :nOfGauss x 1
%}
d = labMesh.d;
X = labMesh.X; % The vector containing the coordinates of the vertices of the labMesh
T = labMesh.T; % The connectivity matrix of the labMesh


nOfGauss = length(locationx); % number of Gauss points in the finite element
depth = zeros(nOfGauss,1);

for g = 1 : nOfGauss
    %[ix,~] = find(X<locationx(g,1));
    [ix,~] = find(X<=locationx(g,1)); 
    labElem = ix(end);
    depth(g) = d(labElem);
end



