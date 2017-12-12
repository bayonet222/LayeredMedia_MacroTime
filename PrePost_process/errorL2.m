%*********************************************************************
%*                        L2 ERROR HDG METHOD                               *
%*********************************************************************
%
%This function calculate the L2 error of the HDG aproximation. 
% L2 error --> SQRT (sum(K in trian) int_K (aprox - exact)^2)
%
%
%***------------------------------------
%*** Inputs: HDG solution : u 
%            Exact solution: uExacta (over each mesh point)
%
%***------------------------------------
% Manuela Bastidas - 2017.

function Error_u = errorL2(X,u,uExacta,geo)


Error_u = 0;

for j=1:geo.nElement
    
    % Coord (x;y) de cada uno de los vertices del tríangulo
    coord = geo.coordinate(geo.element(j,:),:)';
    %     pos = element(j,:); % posiciones de la solución
    pos    = 3*j-2:3*j;
    uAprox = u(pos);
    
    %% Numerical integration
    nCuad = 10;
    [puntos_X1,puntos_Y1,Wx,Wy] = triquad(nCuad,coord');
    test= ([coord;ones(1,3)])\[puntos_X1(:),puntos_Y1(:),ones(nCuad^2,1)]';
    diff = arrayfun(uExacta,repmat(X(1),nCuad^2,1),repmat(X(2),nCuad^2,1),...
        puntos_X1(:),puntos_Y1(:))-test'*uAprox;
    
    
    Error_u = Error_u+ Wx'*reshape(diff.^2,size(Wx,1),size(Wy,1))*Wy;
end

Error_u = sqrt(Error_u);