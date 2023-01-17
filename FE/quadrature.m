function [quadPoints,quadWeights] = quadrature(elemType,numberOfPoints)
%{
---------------------------------------------------------------------------
Input: elemType and numberOfPoints
    elemType      : 0 = line, 1=squares, 2=triangle
    numberOfPoints: number of quadrature points
                   
Output: position and weights of gauss points in the reference element
    quadPoints         : position of gauss points in the reference element dimension 2nu x 2nu
    quadWeights        : weights of gauss points in the reference element with dimension np x 2nu 
    ---------------------------------------------------------------------------
%}
quadPoints = [];
quadWeights = [];
switch (elemType)
    case 0 % 1D lines
        [quadPoints,quadWeights] = lgwt(numberOfPoints,-1,1);
    case 1 % squares
        switch(numberOfPoints)
            case 4 % degree 1
                pos1 = 1/sqrt(3);
                quadPoints = [-pos1 -pos1; pos1 -pos1;
                    pos1  pos1;-pos1  pos1];
                quadWeights = [1 1 1 1];
            case 9 % degree 2 3x3 integration
                pos1 = sqrt(3/5);
                quadPoints = [-pos1 -pos1; 0 -pos1; pos1 -pos1;
                              -pos1     0; 0     0; pos1     0;
                              -pos1  pos1; 0  pos1; pos1  pos1];
                pg1 = 5/9;
                pg2 = 8/9;
                pg3 = pg1;
                quadWeights= [pg1*pg1, pg2*pg1, pg3*pg1, pg1*pg2, pg2*pg2, ...
                    pg3*pg2, pg1*pg3, pg2*pg3, pg3*pg3];
            otherwise
                error('Quadrature for squares with %i points not implemented', numberOfPoints);
        end
    case 2 % triangles
        switch(numberOfPoints)
            case 1 % degree 1
                pos1 = 1/3;
                quadPoints = [pos1  pos1];
                quadWeights = 1/2;
            case 3  % degree 2
                pos1 = 1/2;
                quadPoints = [pos1 pos1; 0 pos1; pos1 0];
                pes1 = 1/6;
                quadWeights = [pes1 pes1 pes1];
            case 4  % degree 3
                quadPoints = [1/3 1/3; 0.6 0.2; 0.2 0.6; 0.2 0.2];
                quadWeights = [-27/96   25/96   25/96   25/96];
            case 7   % degree 5
                a = 0.101286507323456338800987361915123;
                b = 0.470142064105115089770441209513447;
                P1 = 0.0629695902724135762978419727500906;
                P2 = 0.0661970763942530903688246939165759;
                quadPoints = [a a; a 1-2*a; 1-2*a a; b b; b 1-2*b; 1-2*b b; ...
                    1/3 1/3];
                quadWeights = [P1, P1, P1, P2, P2, P2, 0.1125];
            case 12   % degree 6
                a = 0.0630890144915022283403316028708191;
                b = 0.249286745170910421291638553107019;
                c = 0.0531450498448169473532496716313981;
                d = 0.310352451033784405416607733956552;
                P1 = 0.0254224531851034084604684045534344;
                P2 = 0.0583931378631896830126448056927897;
                P3 = 0.0414255378091867875967767282102212;
                quadPoints = [a a; a 1-2*a; 1-2*a a; b b; b 1-2*b; 1-2*b b; ...
                    c d; c 1-c-d; 1-c-d c; 1-c-d d; d 1-c-d; d c];
                quadWeights = [P1, P1, P1, P2, P2, P2, P3, P3, P3, P3, P3, P3];
            case 16  % degree 8
                a = 0.170569307751760206622293501491464;
                b = 0.0505472283170309754584235505965989;
                c = 0.459292588292723156028815514494169;
                d = 0.728492392955404281241000379176061;
                e = 0.263112829634638113421785786284643;
                f = 1/3;
                P1 = 0.0516086852673591251408957751460645;
                P2 = 0.0162292488115990401554629641708902;
                P3 = 0.0475458171336423123969480521942921;
                P4 = 0.0136151570872174971324223450369544;
                P5 = 0.0721578038388935841255455552445323;
                quadPoints = [a a; a 1-2*a; 1-2*a a; b b; b 1-2*b; 1-2*b b;...
                              c c; c 1-2*c; 1-2*c c; d e; d 1-d-e; 1-d-e d; 1-d-e e; ...
                              e 1-d-e; e d; f f];
                quadWeights = [P1*[1 1 1], P2*[1 1 1], P3*[1 1 1], ...
                    P4*ones(1,6), P5];
            otherwise
                error( 'Quadrature for triangles with %i points not implemented', numberOfPoints )
        end
    otherwise
        error( 'Unknown element type: %i \nValid options are 0=1D lines, 1=squares, 2=triangles', elemType )
end


function [x,w] = lgwt(N, a, b)
%
% generates Gauss-Legendre quadratures with N points on the interval [a,b]
%
N  = N-1;  N1 = N+1;   N2 = N+2;

xu = linspace(-1, 1, N1)';

% Initial guess
y = cos((2*(0:N)'+1)*pi/(2*N+2)) + (0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L = zeros(N1, N2);

% Derivative of LGVM
Lp = zeros(N1, N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0 = 2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0)) > eps
    
    L(:,1) = 1;
%     Lp(:,1) = 0;
    
    L(:,2) = y;
%     Lp(:,2) = 1;
    
    for k = 2:N1
        L(:,k+1) = ((2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1))/k;
    end
 
    Lp = (N2)*(L(:,N1)-y.*L(:,N2))./(1-y.^2);
    
    y0 = y;
    y = y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x = (a*(1-y)+b*(1+y))/2;      

% Compute the weights
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

