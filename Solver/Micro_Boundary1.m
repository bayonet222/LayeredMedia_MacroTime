%*********************************************************************
%*                                                                   *
%*                    Micro Boundary constraints                     *
%*                                                                   *
%*********************************************************************
%
% This function is usefull to impose the boundary conditions on the micro
% scale problem.
% Remark: The boundary condition will be Periodic or dirichlet, the user
% can handle it.
% Case Dirichlet: We impose weakle the Dirichlet boundary condition like an
% new data of the liner system.
% Periodic: We impose the peridicity of the soluton at the boundary using
% the selected pairs of edges (Pre-process) and we impose the condition
% like a constraints of the problem, for this reason we use 'lsqlin' to
% solve linear systems with constraints.
%
%***------------------------------------
%***Inputs: PosFlujo, edgesKnum (See Micro_solver)
%           H: Main matrix. J1/2: RHS of each cell problem.
%
%***------------------------------------
% Manuela Bastidas - 2017.

function [SolFlujo1,SolFlujo2] = Micro_Boundary1(posFlujo,...
    edgesKnum,H,J1,J2,Micro_geo,MicroConst)

% global 

nEdgeDir       = Micro_geo.nEdgeDir;
DirichletEdges = Micro_geo.DirichletEdges;
nEdgeTotal     = Micro_geo.nEdgeTotal;

switch MicroConst
    
    case 'Dirichlet'
        
        %% CASE OF THE DIRICHLET CONDITION
        
        IndicadorDirichFlujo = zeros(2*nEdgeTotal,1);
        ValorDirichFlujo     = zeros(2*nEdgeTotal,1);
        
        %################ THIS CAN BE INDEXED
        for edge=1:nEdgeDir
            elemento = DirichletEdges(edge,4);
            pos = posFlujo{elemento,1}...
                (:,edgesKnum{elemento,1}(:,1)==DirichletEdges(edge,1));
            IndicadorDirichFlujo(pos,1) = 1;
            p = Micro_geo.coordinate(DirichletEdges(edge,2:3),:)';
            uD1 = u_D(p(:,1));
            uD3 = u_D(p(:,2));
            pm  = (p(:,1)+p(:,2))./2;
            uD2 = u_D(pm);
            le  = norm(p(:,1)-p(:,2));
            DirDer = le/6*[(uD1+2*uD2);(uD3+2*uD2)];
            ValorDirichFlujo(pos) = (le/3*[1 1/2; 1/2 1])\DirDer;   
        end
        % labels of the position of dirichlet and non-Dirichlet unknowns
        posDir  = find(IndicadorDirichFlujo ==1);
           
        if isempty(posDir)
              posDir = pos;
              % Imposing 0 to the last two degrees of freedom (boundary)
        end
        freePOS = setdiff(1:2*nEdgeTotal,posDir);
        
        % -----------------------------------------------------------
        %            SOLUTION OF DE EDGE PROBLEM
        % -----------------------------------------------------------
        
        SolFlujo1(posDir,1)  = ValorDirichFlujo(posDir);
        SolFlujo1(freePOS,1) = H(freePOS,freePOS)\(J1(freePOS)-...
            H(freePOS,posDir)*ValorDirichFlujo(posDir));
        SolFlujo2(posDir,1)  = ValorDirichFlujo(posDir);
        SolFlujo2(freePOS,1) = H(freePOS,freePOS)\(J2(freePOS)-...
            H(freePOS,posDir)*ValorDirichFlujo(posDir));

    case 'Periodic'
        
        edgePar = Micro_geo.edgePar;
        IndicadorDirichFlujo = eye(2*nEdgeTotal);
        
        Eliminarpos = []; noEliminarpos=[];

        % labels of the related (pairs) of unknowns
        for kk = 1:size(edgePar)
            elemento1 = DirichletEdges((DirichletEdges(:,1)==edgePar(kk,1)),4);
            elemento2 = DirichletEdges((DirichletEdges(:,1)==edgePar(kk,2)),4);
            
            pos1 = posFlujo{elemento1,1}...
                (:,edgesKnum{elemento1,1}(:,1)==edgePar(kk,1));
            pos2 = posFlujo{elemento2,1}...
                (:,edgesKnum{elemento2,1}(:,1)==edgePar(kk,2));
            
            IndicadorDirichFlujo(pos1,pos1) = eye(2);
            IndicadorDirichFlujo(pos2,pos1) = eye(2);
            
            Eliminarpos   = [Eliminarpos pos2'];
            noEliminarpos = [noEliminarpos pos1'];
   
        end
        IndicadorDirichFlujo(:,Eliminarpos)=[];
        
        HH   = (IndicadorDirichFlujo'*H*IndicadorDirichFlujo);
        JJ11 = (IndicadorDirichFlujo'*J1); 
        JJ22 = (IndicadorDirichFlujo'*J2); 
        
        posDir  = [1 2]; %position 1 and 2 imposed to be 0
        freePOS = setdiff(1:length(JJ11),posDir);
        
        % -----------------------------------------------------------
        %            SOLUTION OF DE EDGE PROBLEM
        % -----------------------------------------------------------
        
        Sol1(posDir,1)  = 0;
        Sol1(freePOS,1) = HH(freePOS,freePOS)\JJ11(freePOS);
%         Sol1(freePOS,1) = gmres(HH(freePOS,freePOS),JJ11(freePOS));
        Sol2(posDir,1)  = 0;
        Sol2(freePOS,1) = HH(freePOS,freePOS)\JJ22(freePOS);

        indSolut = setdiff(1:2*nEdgeTotal,full(Eliminarpos));
        SolFlujo1 = zeros(2*nEdgeTotal,1);
        SolFlujo2 = zeros(2*nEdgeTotal,1);
        
        SolFlujo1(indSolut,1)    = Sol1;
        SolFlujo1(Eliminarpos,1) = SolFlujo1(noEliminarpos,1);
        
        SolFlujo2(indSolut,1)    = Sol2;
        SolFlujo2(Eliminarpos,1) = SolFlujo2(noEliminarpos,1);

end

