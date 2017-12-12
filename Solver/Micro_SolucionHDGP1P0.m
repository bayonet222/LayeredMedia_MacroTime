%*********************************************************************
%*                                                                   *
%*                    HDG - Micro CELL problems                      *
%*                                                                   *
%*                  Universiteit Hasselt - CMAT                      *
%*                                                                   *
%*********************************************************************
%
% Solution of the diffusion problem in a Micro-Scale, using linear
% aproximation for the scalar unknown and constant aproximation for the
% vectorial unknown.
%
%***------------------------------------
% Manuela Bastidas - 2017.

function [Micro_Sol] = Micro_SolucionHDGP1P0(Punto_X,tt)

global xi Micro_geo

element      = Micro_geo.element;
coordinate   = Micro_geo.coordinate;

nElement     = Micro_geo.nElement;
% nEdgeNew     = Micro_geo.nEdgeNew;
nEdgeTotal   = Micro_geo.nEdgeTotal;

% nodes2edge   = Micro_geo.nodes2edge;
edge2element = Micro_geo.edge2element;

% NewmanEdges  = Micro_geo.NewmanEdges;

%%
%*********************************************************************
%*                                                                   *
%*                        Inicialization                             *
%*                                                                   *
%*********************************************************************

% First equation of the HDG formulation
A   = cell(nElement,1);
% Second equation of the HDG formulation
C    = cell(nElement,1);
% Matriz RHS
D = cell(nElement,1);
E = cell(nElement,1);
R = cell(nElement,1);

% FLUX Matrix
FF1   = cell(nElement,1);
FF2   = cell(nElement,1);

% Local matrices
M = cell(nElement,1);

X = cell(nElement,1);

Y1 = cell(nElement,1);
Y2 = cell(nElement,1);

% Global matrix
H = sparse(2*nEdgeTotal,2*nEdgeTotal);

K1 = sparse(2*nEdgeTotal,1);
K2 = sparse(2*nEdgeTotal,1);
% Boundary matrix
% G1 = zeros(2*nEdgeTotal,1);

%% Micro Solution

Micro_Sol       = struct();

Micro_Sol.Pres1 = zeros(3*nElement,1);
Micro_Sol.Vel1  = zeros(2*nElement,1);


%% ------------------------------------------------------------------
%                          Position idicators
% -------------------------------------------------------------------

% Vectorial position
posGrad  = Micro_geo.posGrad;
% Scalar position
posPres  = Micro_geo.posPres;
% Numerical Flux position
posFlujo = Micro_geo.posFlujo;
% Sort edges in a triangle
edgesKnum = Micro_geo.edgesKnum;

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

%%
%*********************************************************************
%*                                                                   *
%*                        LOCAL SOLVER 1 - 2                         *
%*                                                                   *
%*********************************************************************

for elemento = 1:nElement
    
    % Coord (x;y) triangle vertex
    coord = coordinate(element(elemento,:),:)';
    
    OrdenposFlujo  = posFlujo{elemento,1}(:);
    normalvector   = Micro_geo.normals{elemento};
    
    nCuad = 4;
    % Integration points - Quadrature
    [X_int,Y_int,Wx,Wy] = triquad(nCuad,coord');
    % Difussion Tensor
    [~,TensorA_inv,RHS_cell] = TensorA_eps(Punto_X,[X_int(:),Y_int(:)],tt);
    
    %***--------------------------------------------------------------
    %                       MATRIZ A
    %***--------------------------------------------------------------
    
    A_aux(1,1) = Wx'*reshape(TensorA_inv(1,1,:),size(Wx,1),size(Wy,1))*Wy;
    A_aux(1,2) = Wx'*reshape(TensorA_inv(1,2,:),size(Wx,1),size(Wy,1))*Wy;
    A_aux(2,1) = Wx'*reshape(TensorA_inv(2,1,:),size(Wx,1),size(Wy,1))*Wy;
    A_aux(2,2) = Wx'*reshape(TensorA_inv(2,2,:),size(Wx,1),size(Wy,1))*Wy;
    A{elemento,1} = A_aux';
    
    %***--------------------------------------------------------------
    %                       MATRIZ F
    %***--------------------------------------------------------------
    
    test= ([coord;ones(1,3)])\[X_int(:),Y_int(:),ones(nCuad^2,1)]';
    F1_aux = repmat(RHS_cell(:,1)',3,1).*test;
    F2_aux = repmat(RHS_cell(:,2)',3,1).*test;
    FF1{elemento,1} = zeros(5,1); FF2{elemento,1} = zeros(5,1);
    for fi = 1:3
        FF1{elemento,1}(2+fi)= Wx'*reshape(F1_aux(fi,:),size(Wx,1),size(Wy,1))*Wy;
        FF2{elemento,1}(2+fi)= Wx'*reshape(F2_aux(fi,:),size(Wx,1),size(Wy,1))*Wy;
    end
    
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
    %                       MATRIZ D
    %***--------------------------------------------------------------
    
    n1 = [normalvector(1,:).*le';normalvector(1,:).*le'];
    n2 = [normalvector(2,:).*le';normalvector(2,:).*le'];
    
    D{elemento,1} = [n1(:)';n2(:)']*1/2;
    
    %***--------------------------------------------------------------
    %                       MATRIZ E
    %***--------------------------------------------------------------
    % ORIENTACIÓN LOCAL DE LAS ARISTAS
    aux = zeros(0);
    aux(element(elemento,:))=1:3;
    aaa = edge2element(edgesKnum{elemento,1},1:2);
    aaa= aux(aaa);
    cambios = find([2 1 1]'-aaa(:,1));
    
    Elocal = zeros(3,6);
    for n=1:3
        if sum(cambios==n)==0
            Eaux1 =  [1/3 1/6; 1/6 1/3]*le(n);
        else
            Eaux1 =  [1/6 1/3; 1/3 1/6]*le(n);
        end
        Elocal([1 2 3]~=n,2*n-1:2*n) = Eaux1;
    end
    E{elemento,1}= xi*Elocal;
    %***--------------------------------------------------------------
    %                       MATRIZ R
    %***--------------------------------------------------------------
    
    R{elemento,1}= xi*[Rlocal(1:2,:)*le(1);
        Rlocal(3:4,:)*le(2); Rlocal(5:6,:)*le(3)];
    
    %%    
    %*****************************************************************
    %                     LOCAL SOLVER 1 - 2
    %*****************************************************************
    
    M{elemento,1}  = [A{elemento,1} zeros(2,3); ...
        zeros(3,2)    xi*C{elemento,1}];
    
    X{elemento,1}  = [D{elemento,1};E{elemento,1}]'*...
        ((M{elemento,1})\[-D{elemento,1};E{elemento,1}])-...
        R{elemento,1};
    
    Y1{elemento,1} = [D{elemento,1};E{elemento,1}]'*...
        ((M{elemento,1})\FF1{elemento,1});
    Y2{elemento,1} = [D{elemento,1};E{elemento,1}]'*...
        ((M{elemento,1})\FF2{elemento,1});
    
    
    H(OrdenposFlujo,OrdenposFlujo) = H(OrdenposFlujo,OrdenposFlujo)+...
        X{elemento,1};
    
    K1(OrdenposFlujo,1) = K1(OrdenposFlujo,1)+Y1{elemento,1};
    K2(OrdenposFlujo,1) = K2(OrdenposFlujo,1)+Y2{elemento,1};
    
end
%%
%*********************************************************************
%*                                                                   *
%*              Solution of the numerical flux's problem             *
%*                                                                   *
%*********************************************************************

[SolFlujo1,SolFlujo2] = Micro_Boundary1(posFlujo,edgesKnum,H,-K1,-K2);

%%
%*********************************************************************
%*                                                                   *
%*                   Solution of each LOCAL SOLVER                   *
%*                                                                   *
%*********************************************************************

for elemento = 1:nElement
    posflujoLocal = posFlujo{elemento,1}(:);
    
    RHSlocal1 = [-D{elemento,1};E{elemento,1}]*...
        SolFlujo1(posflujoLocal,1)+FF1{elemento,1};
    RHSlocal2 = [-D{elemento,1};E{elemento,1}]*...
        SolFlujo2(posflujoLocal,1)+FF2{elemento,1};
    
    % Deafault Matlab Linear system solver
    solLocal1 = M{elemento,1}\RHSlocal1;
    solLocal2 = M{elemento,1}\RHSlocal2;
    
    [~,invTensorA,~] = TensorA_eps(Punto_X,sum(coordinate(element(elemento,:),:))/3,tt);
    
    % Saving results
    Micro_Sol.Pres1(posPres{elemento,1},1) = solLocal1(3:end);
    Micro_Sol.Pres2(posPres{elemento,1},1) = solLocal2(3:end);
    
    Micro_Sol.Vel1(posGrad{elemento,1},1)  = -invTensorA*solLocal1(1:2);
    Micro_Sol.Vel2(posGrad{elemento,1},1)  = -invTensorA*solLocal2(1:2);

end


