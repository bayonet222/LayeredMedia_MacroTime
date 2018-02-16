%*********************************************************************
%*                                                                   *
%*                    HDG - MACRO UPSCALED PROBLEM                   *
%*                                                                   *
%*                  Universiteit Hasselt - CMAT                      *
%*                                                                   *
%*********************************************************************
%
% Solution of the upscaled diffusion problem (Macro), using linear
% aproximation for the scalar unknown and constant aproximation for the
% vectorial unknown.
%
%***------------------------------------
%***Inputs: A_efective -> 0 indicates that the diffusion tensor does
%           not depend of macro scale.
%
%***------------------------------------
% Manuela Bastidas - 2017.


function [Macro_tSol] = Macro_SolucionHDGP1P1(Pre_ant,matAux,A_Efective,...
    Pre_Aefect,t_pos,Macro_geo,Adepend,Time,Micro_geo)

% global  Macro_geo  xi Adepend Time
 xi = 1;
element        = Macro_geo.element;
coordinate     = Macro_geo.coordinate;

nElement       = Macro_geo.nElement;
nEdgeDir       = Macro_geo.nEdgeDir;
nEdgeTotal     = Macro_geo.nEdgeTotal;

% nodes2edge     = Macro_geo.nodes2edge;
% edge2element   = Macro_geo.edge2element;

DirichletEdges = Macro_geo.DirichletEdges;

%%
%*********************************************************************
%*                                                                   *
%*                        Inicialization                             *
%*                                                                   *
%*********************************************************************

% First equation of the HDG formulation
A   = cell(nElement,1);
B   = matAux.B;
% Second equation of the HDG formulation
C    = matAux.C;

% Matriz RHS
D = matAux.D;

E = matAux.E;

R = matAux.R;
T = matAux.T;

% FLUX Matrix
FF1   = cell(nElement,1);

% Local matrices
M = cell(nElement,1);

X = cell(nElement,1);

Y1 = cell(nElement,1);

% Global matrices
H  = sparse(2*nEdgeTotal,2*nEdgeTotal);
K1 = sparse(2*nEdgeTotal,1);

% Auxiliar Bpundary matrices
IndicadorDirichFlujo = zeros(2*nEdgeTotal,1);
ValorDirichFlujo     = zeros(2*nEdgeTotal,1);


%% Micro Solution

Macro_tSol       = struct();

Macro_tSol.Pres  = zeros(3*nElement,1);
Macro_tSol.Vel   = zeros(6*nElement,1);
tt = Time.time_vec(t_pos);

%% ------------------------------------------------------------------
%                          Position idicators
% -------------------------------------------------------------------

% Vectorial position
posGrad  = Macro_geo.posGrad;
% Scalar position
posPres  = Macro_geo.posPres;
% Numerical Flux position
posFlujo = Macro_geo.posFlujo;
% Sort edges in a triangle
edgesKnum = Macro_geo.edgesKnum;


%% ------------------------------------------------------------------
%                          Auxiliar matrices idicators
% -------------------------------------------------------------------

Alocal = 1/12*[1 1/2 1/2 ;
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
    OrdenposFlujo  = posFlujo{elemento,1}(:);
    %     normalvector   = Macro_geo.normals{elemento};
    
    %***--------------------------------------------------------------
    %                       MATRIZ A
    %***--------------------------------------------------------------
    nCuad = 2;
    [puntos_X1,puntos_Y1,Wx,Wy] = triquad(nCuad,coord'); %N^2 puntos
    % coefK: [a_j ,b_j, c_j]
    coefK = inv([coord;ones(1,3)]);
    test= coefK*[puntos_X1(:),puntos_Y1(:),ones(nCuad^2,1)]';
    
    
    %% Micro cell problem for each Macro point
    
    %     [A{elemento,1},A_Efective_mean]= Macro_matA_localsolver(puntos_X1,puntos_Y1,Wx,Wy,...
    %     test,1);
    if strcmp(Adepend,'Macro')
        
        if Time.MicroTdepend == 1 || t_pos == 2 %first time step

%             fprintf('\n Time %i/%i - Element %i/%i \n',...
%                 t_pos,Time.tnSteps,elemento,nElement)
            
            MacroP = sum(coordinate(element(elemento,:),:))/3;
            [A{elemento,1},A_Efective_mean]= Macro_matA_localsolver2(Wx,Wy,...
                test,1,MacroP,tt,Micro_geo);
            
            Macro_tSol.A_Efective(:,:,elemento)= A_Efective_mean;
            
        else
            
            areaK = Macro_geo.areaK{elemento,1};
            A{elemento,1} = areaK.*kron(inv(Pre_Aefect(:,:,elemento))',Alocal)*2;
            Macro_tSol.A_Efective(:,:,elemento)= Pre_Aefect(:,:,elemento);
        end
    else
        areaK = Macro_geo.areaK{elemento,1};
        A{elemento,1} = areaK.*kron(inv(A_Efective)',Alocal)*2;
    end
    %     Macro_Sol.Microerror1(elemento) = max(microerror(1,:));
    %     Macro_Sol.Microerror2(elemento) = max(microerror(2,:));
    %
    
    %***--------------------------------------------------------------
    %                       MATRIZ F
    %***--------------------------------------------------------------
    
    Faux = arrayfun(@f,puntos_X1(:),puntos_Y1(:),repmat(tt,nCuad^2,1));
    
    Faux = repmat(Faux,1,3).*test';
    FF1{elemento,1} = zeros(9,1);
    for fi = 1:3
        FF1{elemento,1}(6+fi)= Wx'*reshape(Faux(:,fi),size(Wx,1),size(Wy,1))*Wy;
    end
    
    
    %*****************************************************************
    %                     LOCAL SOLVER
    %*****************************************************************
    
    M{elemento,1}  = [A{elemento,1} -B{elemento,1}; ...
        Time.dt*B{elemento,1}'   T{elemento,1} + Time.dt*xi*C{elemento,1}];
    
    X{elemento,1}  = [D{elemento,1};E{elemento,1}]'*...
        ((M{elemento,1})\[-D{elemento,1};Time.dt*E{elemento,1}])-...
        xi*R{elemento,1};
    
    Y1{elemento,1} = [D{elemento,1};E{elemento,1}]'*...
        ((M{elemento,1})\(Time.dt*FF1{elemento,1}));
    Y1{elemento,1} = Y1{elemento,1} + ...
        [D{elemento,1};E{elemento,1}]'*...
        ((M{elemento,1})\[zeros(6,1);(T{elemento,1}*Pre_ant(posPres{elemento,1},1))]);
    
    H(OrdenposFlujo,OrdenposFlujo) = H(OrdenposFlujo,OrdenposFlujo)+...
        X{elemento,1};
    
    K1(OrdenposFlujo,1) = K1(OrdenposFlujo,1)+Y1{elemento,1};
end

%%
%*********************************************************************
%*                                                                   *
%*                        DIRICHLET BOUNDARY                           *
%*                                                                   *
%*********************************************************************

for edge=1:nEdgeDir
    
    elemento = DirichletEdges(edge,4);
    
    pos = posFlujo{elemento,1}(:,edgesKnum{elemento,1}(:,1)==DirichletEdges(edge,1));
    IndicadorDirichFlujo(pos,1) = 1;
    
    %     p = Macro_geo.coordinate(DirichletEdges(edge,2:3),:)';
    %     uD1 = u_D(p(:,1));
    %     uD3 = u_D(p(:,2));
    %     pm  = (p(:,1)+p(:,2))./2;
    %     uD2 = u_D(pm);
    %     le = norm(p(:,1)-p(:,2));
    %     DirDer = le/6*[(uD1+2*uD2);(uD3+2*uD2)];
    %
    %     ValorDirichFlujo(pos) = (le/3*[1 1/2; 1/2 1])\DirDer;
end

%%
%*********************************************************************
%*                                                                   *
%*              Solution of the numerical flux's problem             *
%*                                                                   *
%*********************************************************************

% Imposing dirichlet condition - Tmp matrices
posDir  = find(IndicadorDirichFlujo ==1);
freePOS = setdiff(1:2*nEdgeTotal,posDir);

SolFlujo1(posDir,1)  = ValorDirichFlujo(posDir);
SolFlujo1(freePOS,1) = H(freePOS,freePOS)\(-K1(freePOS)-...
    H(freePOS,posDir)*ValorDirichFlujo(posDir));


%%
%*********************************************************************
%*                                                                   *
%*                   Solution of each Macro LOCAL SOLVER             *
%*                                                                   *
%*********************************************************************

for elemento = 1:nElement
    posflujoLocal = posFlujo{elemento,1}(:);
    
    RHSlocal1 = [-D{elemento,1};Time.dt*E{elemento,1}]*...
        SolFlujo1(posflujoLocal,1)+Time.dt*FF1{elemento,1}+...
        [zeros(6,1);(T{elemento,1}*Pre_ant(posPres{elemento,1},1))];
    
    % Deafault Matlab Linear system solver
    solLocal1 = M{elemento,1}\RHSlocal1;
    
    % Saving results
    Macro_tSol.Pres(posPres{elemento,1},1)  = solLocal1(7:end);
    
    if strcmp(Adepend,'Macro')
        if t_pos==2
            IA_Efective = inv(Macro_tSol.A_Efective(:,:,elemento));
        else
            IA_Efective = inv(Pre_Aefect(:,:,elemento));
        end
    else
        IA_Efective = inv( A_Efective);
    end
    
    %posprocess velocity
    Solvel_aux = zeros(6,1);
    for kk=1:3
        Solvel_aux([kk,3+kk])= -IA_Efective*solLocal1([kk,3+kk]);
    end
    Macro_tSol.Vel(posGrad{elemento,1},1) = Solvel_aux;
end


