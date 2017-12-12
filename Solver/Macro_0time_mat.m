%*********************************************************************
%*                                                                   *
%*            HDG - MACRO UPSCALED PROBLEM - Aux Matrices            *
%*                                                                   *
%*********************************************************************
%
% Solution of the upscaled diffusion problem (Macro), using linear
% aproximation for the scalar unknown and constant aproximation for the
% vectorial unknown.
% Auxiliar code to compute the matrices that does not depend of the time
% step
%
%***------------------------------------
% Manuela Bastidas - 2017.


function [Macro_0time_matAux] = Macro_0time_mat

global  Macro_geo xi

element        = Macro_geo.element;
coordinate     = Macro_geo.coordinate;

nElement       = Macro_geo.nElement;

% nodes2edge     = Macro_geo.nodes2edge;
edge2element   = Macro_geo.edge2element;

%%
%*********************************************************************
%*                                                                   *
%*                        Inicialization                             *
%*                                                                   *
%*********************************************************************


% Second equation of the HDG formulation
B    = cell(nElement,1);
C    = cell(nElement,1);

% Matriz RHS
D = cell(nElement,1);

E = cell(nElement,1);

R = cell(nElement,1);
T = cell(nElement,1);

%% ------------------------------------------------------------------
%                          Position idicators
% -------------------------------------------------------------------

% Vectorial position
% posGrad  = Macro_geo.posGrad;
% % Scalar position
% posPres  = Macro_geo.posPres;
% % Numerical Flux position
% posFlujo = Macro_geo.posFlujo;
% Sort edges in a triangle
edgesKnum = Macro_geo.edgesKnum;



%% ------------------------------------------------------------------
%                          Auxiliar matrices idicators
% -------------------------------------------------------------------

% Local R matrix in ref triangle.
Rlocal = 1/3*[1 1/2 0 0 0 0 ;
    1/2 1 0 0 0 0;
    0 0 1 1/2 0 0;
    0 0 1/2 1 0 0;
    0 0 0 0 1 1/2;
    0 0 0 0 1/2 1];

Tlocal = 1/12*[1 1/2 1/2 ;
    1/2 1  1/2  ;
    1/2 1/2 1] ;

%%
%*********************************************************************
%*                                                                   *
%*                        LOCAL SOLVER (MACRO)                         *
%*                                                                   *
%*********************************************************************

for elemento = 1:nElement
    
    % Coord (x;y) triangle vertex
    coord = coordinate(element(elemento,:),:)';
    normalvector   = Macro_geo.normals{elemento};
    
    %***--------------------------------------------------------------
    %                       MATRIZ B
    %***--------------------------------------------------------------
    nCuad = 4;
    [puntos_X1,puntos_Y1,Wx,Wy] = triquad(nCuad,coord'); %N^2 puntos
    % coefK: [a_j ,b_j, c_j]
    coefK = inv([coord;ones(1,3)]);
    test= coefK*[puntos_X1(:),puntos_Y1(:),ones(nCuad^2,1)]';
    testMat = zeros(6,3);
    
    for kk = 1:3
        testMat(:,kk) = ones(6,1).*(Wx'*reshape(test(kk,:),size(Wx,1),size(Wy,1))*Wy);
    end
    B{elemento,1} = testMat.*repmat([coefK(:,1);coefK(:,2)],1,3);
    
     
    %***--------------------------------------------------------------
    %                       MATRIZ C
    %***--------------------------------------------------------------
    % coord egdes : [1:2] inicial , [3:4] final
    p      = [coord(:,[2 3 1]);coord(:,[3 1 2])];
    % Length edge
    le     = [norm(p(1:2,1)-p(3:4,1));norm(p(1:2,2)-p(3:4,2));...
        norm(p(1:2,3)-p(3:4,3))];
    
    C{elemento,1} = 1/3*[le(2)+le(3) le(3)/2 le(2)/2;...
        le(3)/2 le(1)+le(3) le(1)/2;...
        le(2)/2 le(1)/2 le(1)+le(2)];
    
    %***--------------------------------------------------------------
    %                       MATRIZ E
    %***--------------------------------------------------------------
    
    % ORIENTACIÓN LOCAL DE LAS ARISTAS
    aux = zeros(0);
    aux(element(elemento,:))=1:3;
    aaa = aux(edge2element(edgesKnum{elemento,1},1:2));
    cambios = find([2 1 1]'-aaa(:,1));
    
    for n=1:3
        if sum(cambios==n)==0
            Eaux1 =  [1/3 1/6; 1/6 1/3]*le(n);
        else
            Eaux1 =  [1/6 1/3; 1/3 1/6]*le(n);
        end
        E{elemento,1}([1 2 3]~=n,2*n-1:2*n) = xi*Eaux1;
    end
    
    n1 = [normalvector(1,:);normalvector(1,:)];
    Dlocal(1,:) = n1(:)';
    DlocalNorm(1:3,:) = repmat(Dlocal(1,:),3,1);
    n2 = [normalvector(2,:);normalvector(2,:)];
    Dlocal(2,:) = n2(:)';
    DlocalNorm(4:6,:) = repmat(Dlocal(2,:),3,1);
    
    D{elemento,1} = [E{elemento,1};E{elemento,1}].*DlocalNorm;
    
    
    %***--------------------------------------------------------------
    %                       MATRIZ R
    %***--------------------------------------------------------------
    
    R{elemento,1}= [Rlocal(1:2,:)*le(1);
        Rlocal(3:4,:)*le(2); Rlocal(5:6,:)*le(3)];
    
   %***--------------------------------------------------------------
   %                       MATRIZ T - Diffusion Problem
   %***--------------------------------------------------------------
    
    T{elemento,1} = Macro_geo.areaK{elemento,1}*Tlocal*2;
end

Macro_0time_matAux.B = B;
Macro_0time_matAux.T = T;
Macro_0time_matAux.C = C;
Macro_0time_matAux.D = D;
Macro_0time_matAux.E = E;
Macro_0time_matAux.R = R;


