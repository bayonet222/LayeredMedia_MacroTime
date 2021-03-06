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

function Error_u = errorLinf_Macro(t_pos,u,uExacta,geo)

global Time
tt = Time.time_vec(t_pos);

Error_u = 0;

for j=1:geo.nElement
    
    % Coord (x;y) de cada uno de los vertices del tr�angulo
    coord = geo.coordinate(geo.element(j,:),:)';
    %     pos = element(j,:); % posiciones de la soluci�n
    pos    = 3*j-2:3*j;
    uAprox = u(pos);
    
    %% Numerical integration
    nCuad = 10;
    [puntos_X1,puntos_Y1,~,~] = triquad(nCuad,coord');
    test= ([coord;ones(1,3)])\[puntos_X1(:),puntos_Y1(:),ones(nCuad^2,1)]';
    diff = arrayfun(uExacta,puntos_X1(:),puntos_Y1(:),...
        repmat(tt,nCuad^2,1))-test'*uAprox;
    
    Error_u = max(Error_u,max(diff));
end

% Error_u = sqrt(Error_u);