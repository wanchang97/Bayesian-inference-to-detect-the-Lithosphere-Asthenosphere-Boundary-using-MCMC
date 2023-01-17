function depth = LABdepth2(locationx,labMesh)
%{
To find the LABdepth of all gauss points in the Finite element
Here the code is for the (Nearest Neighbor)Voronoin interpolation,discrete
For future we can have (linear)Delaunay interpolation or natural neighbor
(Voronois +  weighted average)
Input: 
    d         ------ LAB depth discretization |dim : np x 1
    locationx ------ the given location | dim : (nsd-1) x nOfGauss 
    labMesh   .C    ------ LAB depth discretization |dim : np x (nsd-1)
              .type ------ LAB depth discretization |dim : '1D' or '2D'
              .T    ------ LAB depth discretization |dim : np x ()
              .X    ------ LAB depth discretization |dim : (np+1) x (nsd-1)

Output: 
    depth     ------ the LAB depth of all the mid points
                                         |dim :nOfGauss x 1
%}
d = labMesh.d;
X = labMesh.X(2:end); % The vector containing the coordinates of the vertices of the labMesh
%T = labMesh.T; % The connectivity matrix of the labMesh
nOfelements= size(locationx,2); % number of Gauss points in the finite element
depth = zeros(nOfelements,1);
switch labMesh.type
    case '1D'
        for e = 1 : nOfelements
            [ix,~] = find(X>locationx(:,e));
            labElem = ix(1);
            depth(e) = d(labElem);
        end
    case '2D'
        
    otherwise
        error('LABmesh tydpe: 1D or 2D.')
end


% d = reshape(d,1,[]);
% % whichcolumn  = fix(locationx ./param.columnWidth)+1; 
% % depth = m(:,whichcolumn);
% % for 1D discretization of LAB: 
% X = labMesh.X;
% %T = labMesh.T;
% n = length(locationx);
% %np = size(T,1); % number of elements
% 
% depth = zeros(n,1);
% for g = 1:n
%     [~,ix] = find(X<locationx(g));
%     %index = find(X>=locationx,'first');
%     labElem = ix(end);
%     depth(g) = d(labElem);   
% end
% end



