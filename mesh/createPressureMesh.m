function [XP,TP] = createPressureMesh(elemType,nenP,X,T,nx,ny)
%{
---------------------------------------------------------------------------
The file could provide us the coordinate matrix and the connectivity matrix
for pressure
---------------------------------------------------------------------------
Input: 
    elemType  : 1=squares Q, 2=triangles P
    nenP      : number of pressure nodes in the element
    X         : The coordinate matrix for velocity contains the information of global nodes
                numbering and the coordinates with dimension nu x 2
    T         : The Connectivity matrix for velocity contains the information of the element
                and the global nodes numbering written in anticlock wise
                with dimension: number of elements x number of nodes per
                element
    nx        : number of elements in x direction
    ny        : number of elements in y direction

                   
Output: 
    XP        : The coordinate matrix for pressure XP stacks all the coordinates of the nodes in a vector in 
                a sequence of "the global nodes numbering" with dimension np x 2
    TP        : The Connectivity matrix for pressure TP contains the element
                and "the global nodes numbering" written in anticlock wise
                with dimension: number of elements x number of nodes per
                element
---------------------------------------------------------------------------           
%}
[numel,nen] = size(T); 

switch(elemType)
   case 1 % squares
      switch(nenP)
         case 1 % P0
            XP = zeros(numel,2); TP = zeros(numel,nenP);
            for ielem=1:numel
               XP(ielem,:) = mean(X(T(ielem,:),:));
               TP(ielem)   = ielem;
            end
            
         case 4
            if nen == 4  % Q1Q1
               XP = X; TP = T;
            elseif  nen==9 % Q2Q1
               % regular structures meshes assumed!
               npx = nx/2+1; 
               npy = ny/2+1;
               XP = zeros(npx*npy,2); 
               TP = zeros(numel,nenP);
               for irow = 1:npy
                  irowodd = 2*(irow-1)+1;
                  XP((irow-1)*npx+1:irow*npx,:) = X((irowodd-1)*(nx+1)+1:2:irowodd*(nx+1),:);
               end
               for i = 1:ny/2
                  for j = 1:nx/2
                     ielem = (i-1)*nx/2+j;
                     inode = (i-1)*(npx)+j;
                     TP(ielem,:) = [inode inode+1 inode+(npx+1) inode+(npx)];
                  end
               end
            end
            
         case 9
            if nen==9 % Q2Q2
               XP = X;
               TP = T;
            end
      end
      
   case 2 % Triangles
      switch(nenP)
         case 1 % P0
            XP = zeros(numel,2); 
            TP = zeros(numel,nenP);
            for ielem = 1:numel
               XP(ielem,:) = mean(X(T(ielem,:),:));
               TP(ielem) = ielem;
            end
            
         case 3
            if nen == 3 % P1P1
               XP = X;
               TP = T;
               
            elseif nen == 4 % MINI (P1+P1)
               XP = X(1:(nx+1)*(ny+1),:);
               TP = zeros(nx*ny,3); 
               TP = T(:,1:3);
               
            elseif nen == 6 % P2P1
               npx = nx/2+1; 
               npy = ny/2+1;
               XP = zeros(npx*npy,2);
               TP = zeros(nx*ny/4,3);
               for irow = 1:npy
                  irowodd = 2*(irow-1)+1;
                  XP((irow-1)*npx+1:irow*npx,:) = X((irowodd-1)*(nx+1)+1:2:irowodd*(nx+1),:);
               end
               for i = 1:ny/2
                  for j = 1:nx/2
                     ielem = 2*((i-1)*(nx/2)+j)-1;
                     inode = (i-1)*(npx)+j;
                     TP(ielem,:) = [inode   inode+1   inode+(npx)];
                     TP(ielem+1,:) = [inode+1   inode+1+npx   inode+npx];
                  end
               end
               TP(1,:) = [1  npx+2   npx+1];
               TP(2,:) = [1     2    npx+2];
               aux = size(TP,1);
               TP(aux,:)   = [npx*(npy-1)-1   npx*(npy-1) npx*npy];
               TP(aux-1,:) = [npx*(npy-1)-1   npx*npy     npx*npy-1];
            end
         case 6
            
            if nen == 6  % P2P2
               XP = X;
               TP = T;
            end
      end
end

