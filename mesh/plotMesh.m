function h = plotMesh(X,T,varargin)
%{
---------------------------------------------------------------------------
Plot the mesh given the coordinate matrix X and the connectivity matrix T
---------------------------------------------------------------------------
Input: 
    X         : The coordinate matrix with information of global
                 coordinates and global numbering of the nodes 121x2
    T         : The coordinate matrix with information of global and 
                local nodes numebring and the element numbering 25x9
    varargin  :               
        h = plotMesh(...,'plotNodes',1)
        h = plotMesh(...,'lineSpecNodes','.k')
        h = plotMesh(..., 'labelNodes', 1 )
        h = plotMesh(..., 'labelNodesFontSize', 8)
        h = plotMesh(..., 'labelNodesColor', 'k')
        h = plotMesh(..., 'plotElements', 1)
        h = plotMesh(..., 'lineSpecElements', '-b')
        h = plotMesh(..., 'lineThicknessElements', 0.5)
        h = plotMesh(..., 'labelElements', 1)
        h = plotMesh(..., 'labelElementsFontSize', 8)
        h = plotMesh(..., 'labelElementsColor', 'b')
        h = plotMesh(..., 'elementType', e)
            e = 0 line
            e = 1 square
            e = 2 triangle
Output: 
    h:        : figure handle to create objects
---------------------------------------------------------------------------           
%}

%X = mesh.XP;
%T = mesh.TP;

plotNodes = getValueFromVarargin('plotNodes',varargin,0);
lineSpecNodes = getValueFromVarargin('lineSpecNodes',varargin,'.k');
labelNodes = getValueFromVarargin('labelNodes',varargin,0);
labelNodesFontSize = getValueFromVarargin('labelNodesFontSize',varargin,8);
labelNodesColor = getValueFromVarargin('labelNodesColor',varargin,'k');
plotElements = getValueFromVarargin('plotElements',varargin,1);
lineSpecElements = getValueFromVarargin('lineSpecElements',varargin,'-b');
lineThicknessElements = getValueFromVarargin('lineThicknessElements',varargin,0.5);
labelElements = getValueFromVarargin('labelElements',varargin,0);
labelElementsFontSize = getValueFromVarargin('labelElementsFontSize',varargin,8);
labelElementsColor = getValueFromVarargin('labelElementsColor',varargin,'b');
elemType = getValueFromVarargin('elementType',varargin);
%
dims = size( X, 2 );
hnodes = [];
helems = [];
if isempty(elemType)
   switch(size(T,2))
      case 2
         elemType = 0; % linear 1D
      case 3
         elemType = 2; % linear triangle
      case 4
         maxnode = max(T(:));
         if maxnode > size(X,1)
            elemType = 2; % mini triangular elem
         else
            elemType = 1; % linear square
         end
      otherwise
         error('element not implemented')
   end
end

switch(elemType)
   case 2
      switch( size(T,2) )
         case 1
            ix = 1;
         case {3,4,6}
            ix = [1 2 3 1];
      end
   case 1
      switch( size(T,2) )
         case 1
            ix = 1;
         otherwise
            ix = [1 2 3 4];
      end
end
hold on
if plotNodes && ~labelNodes
   if dims == 2
      hnodes = plot( X(:,1), X(:,2), lineSpecNodes, ...
         'Tag', 'nodeLineTag' );
   else
      hnodes = plot3( X(:,1), X(:,2), X(:,3),lineSpecNodes, ...
         'Tag', 'nodeLineTag' );
   end
end
if labelNodes
   hnodes = zeros( 1, length( X ) );
   for I = 1:length( X )
      if dims == 2
         hnodes(I) = text( X(I,1), X(I,2), int2str( I ), ...
            'FontSize', labelNodesFontSize, ...
            'Color', labelNodesColor, ...
            'Tag', 'nodeNumberTag' );
      else
         hnodes(I) = text( X(I,1), X(I,2), X(I,3), int2str( I ), ...
            'FontSize', labelNodesFontSize, ...
            'Color', labelNodesColor, ...
            'Tag', 'nodeNumberTag' );
      end
   end
end
if plotElements && ~isempty( T )
   if dims==2
      helems = zeros( 1, length( T ) );
      for j = 1:size( T, 1 )
         helems(j) = plot( X(T(j,ix),1), X(T(j,ix),2), ...
            lineSpecElements, ...
            'lineWidth', lineThicknessElements, ...
            'Tag', 'elementLineTag' );
         if labelElements
            prom = mean( X( T(j,1:3),:) );
            helems(j) = text( prom(1), prom(2), int2str( j ), ...
               'FontSize', labelElementsFontSize, ...
               'Color', labelElementsColor, ...
               'Tag', 'elementNumberTag' );
         end
      end
   else
      helems = [];
      warning('common:plot:plotMesh', 'we do not plot 3D elements')
   end
end
h = [hnodes helems];
