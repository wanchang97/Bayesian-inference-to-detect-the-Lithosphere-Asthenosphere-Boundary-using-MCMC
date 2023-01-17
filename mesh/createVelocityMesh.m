function [X,T] = createVelocityMesh(elemType,nen,x1,x2,y1,y2,nx,ny)
%{
---------------------------------------------------------------------------
The file could provide us the coordinate matrix and the connectivity matrix for velocity
---------------------------------------------------------------------------
Input: 
    elemType  : 1=squares Q, 2=triangles P
    nen       : number of velocity nodes in the element
    x1        : domain size xmin
    x2        : domain size xmax
    y1        : domain size ymin
    y2        : domain size ymax
    nx        : number of elements in x direction
    ny        : number of elements in y direction

                   
Output: 
    X         : The coordinate matrix for velocity X stacks all the coordinates of the nodes in a vector in 
                a sequence of "the global nodes numbering" with dimension nu x 2
    T         : The Connectivity matrix for velocity T contains the element
                and "the global nodes numbering" written in anticlock wise
                with dimension: number of elements x number of nodes per
                element
---------------------------------------------------------------------------           
%}
% number of nodes in each direction
npx = nx+1; 
npy = ny+1;

X = zeros((npx)*(npy),2);
xs = linspace(x1,x2,npx)'; 

% nodal coords
yys = linspace(y1,y2,npy);
for i=1:npy
 ys = yys(i)*ones(npx,1);
 posi = (i-1)*(npx)+1:i*(npx);
 X(posi,:)=[xs,ys];
end

switch(elemType)
   case 1 % squares
      switch(nen)
         case 4 % Q1
            T = zeros(nx*ny,4);
            for i = 1:ny
               for j = 1:nx
                  ielem = (i-1)*nx+j;
                  inode = (i-1)*(npx)+j;
                  T(ielem,:) = [inode inode+1 inode+(nx+2) inode+(npx)];
               end
            end
         case 9 % Q2
            if (nx-2*floor(nx/2) ~=0) && (ny-2*floor(ny/2) ~=0)
               error('Number of nodes in X or Y direction is not odd')
            end
            for i = 1:ny/2
               for j = 1:nx/2
                  ielem=(i-1)*nx/2+j;
                  inode=(i-1)*2*(npx)+2*(j-1)+1;
                  T(ielem,:)=[inode inode+2 inode+2*(npx)+2 inode+2*(npx) inode+1 inode+2+(npx) inode+2*(npx)+1 inode+(npx) inode+(npx)+1];
               end
            end
      end
      
   case 2 % Triangles
      switch(nen)
         case 3 % P1
            T = zeros(nx*ny,3);
            for i = 1:ny
               for j = 1:nx
                  ielem = 2*((i-1)*nx+j)-1;
                  inode = (i-1)*(npx)+j;
                  T(ielem,:) = [inode   inode+1   inode+(npx)];
                  T(ielem+1,:) = [inode+1   inode+1+npx   inode+npx];
               end
            end
            % Change the elements of the SW and NE corners to avoid
            % prescribing all dofs in the element
            T(1,:) = [1  npx+2   npx+1];
            T(2,:) = [1    2     npx+2];
            aux = size(T,1);
            T(aux,:) = [npx*ny-1  npx*npy  npx*npy-1];
            T(aux-1,:) = [npx*ny-1  npx*ny  npx*npy];
            
         case 4 % P1+ (mini) the interior node is NOT generated
            T = zeros(nx*ny,4);
            for i = 1:ny
               for j = 1:nx
                  ielem = 2*((i-1)*nx+j)-1;
                  inode = (i-1)*(npx)+j;
                  n_ad = npx*npy + 2*((i-1)*nx+j)-1;
                  T(ielem,:) = [inode  inode+1  inode+npx  n_ad];
                  T(ielem+1,:) = [inode+1  inode+1+npx  inode+npx  n_ad+1];
               end
            end
            % Change the elements of the SW and NE corners to avoid
            % prescribing all dofs in the element
            aux = size(T,1);
            T(1,:) = [1  npx+2   npx+1  npx*npy+1];
            T(2,:) = [1      2   npx+2  npx*npy+2];
            T(aux-1,:) = [npx*ny-1  npx*npy  npx*npy-1 npx*npy+2*nx*ny-1];
            T(aux,:)   = [npx*ny-1  npx*ny   npx*npy   npx*npy+2*nx*ny];
            
         case 6 % P2
            for i = 1:ny/2
               for j = 1:nx/2
                  ielem = 2*((i-1)*nx/2+j) - 1;
                  inode = (i-1)*2*(npx)+2*(j-1)+1;
                  T(ielem,:) = [inode   inode+2   inode+2*npx   inode+1    inode+1+npx   inode+npx];
                  T(ielem+1,:) = [inode+2    inode+2+2*npx   inode+2*npx   inode+2+npx   inode+1+2*npx   inode+1+npx];
               end
            end
            % Change the elements of the SW and NE corners to avoid
            % prescribing all dofs in the element
            T(1,:) = [1   2*npx+3    2*npx+1   npx+2   2*npx+2   npx+1];
            T(2,:) = [1      3       2*npx+3     2      npx+3    npx+2];
            aux = size(T,1);
            T(aux-1,:) = [npx*(ny-1)-2   npx*(ny-1)    npx*npy    npx*(ny-1)-1    npx*ny     npx*ny-1];
            T(aux,:)   = [npx*(ny-1)-2    npx*npy     npx*npy-2     npx*ny-1     npx*npy-1   npx*ny-2 ];
      end
end

