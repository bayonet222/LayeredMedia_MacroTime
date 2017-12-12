%*********************************************************************
%*                                                                   *
%*                      HDG - Upscaling                              *
%*                                                                   *
%*                 Manuela Bastidas  - 2017                          *
%*               Universiteit Hasselt - CMAT                         *
%*                                                                   *
%*********************************************************************
%
% Two-dimensional  HDG method for elliptic problems.
%
% HDG_UP solves the equation
%     - div(a*grad(u)) = f in Omega
%                    u = 0 on the Dirichlet boundary
%               d/dn u = 0   on the Neumann boundary
% with tensor a(x) on a geometry described by trian.
%
% See folder Inputs
%
%***------------------------------------
% Manuela Bastidas  - 2017.

%%
%*********************************************************************
%*                                                                   *
%*                      Simulation Parameters                        *
%*                                                                   *
%*********************************************************************

% Preprocess time start
% tic;
%
clear
close all
clc


global Adepend MicroConst
% Adepend -> Dependency of the tensor A
% Adepend = Micro : A(X,Y) = A(.,Y) (Micro dependency)
% Adepend = Macro : A(X,Y) = A(X,Y) (Macro-Micro dependency)

% MicroConst = Dirichlet %-> Dirichlet boundary conditions - MicroScale
% MicroConst = Periodic  %-> Periodic boundary conditions - MicroScale

% Adepend    = 'Macro';
% MicroConst = 'Periodic';

% % Name of the Gmsh files
% % 1 -> MACRO MESH
% % 2 -> MICRO MESH : geo_MicroRef REGULAR MESH OF [0,1]^2

% meshfolder  = '\Mesh\UnitSquare\';
% file1 = '\UnitSquare_20.msh';
% file2 = '\UnitSquare_20.msh';

% % Boundary labels
% Macro Boundary labels
% indGamma_N1  = NaN;      % Newmann
% indGamma_D1  = 7 ;       % Dirichlet

% Micro Boundary Labels
% indGamma_N2  = NaN;      % Newmann   -> Grain Boundary
% indGamma_D2  = 7;        % Dirichlet -> External micro boundary

% --------------------------------------------------------------------

% Files to pre - process (read Mesh, create mesh object and auxiliars)
addpath([cd,'\PrePost_process'])
addpath([cd,'\Inputs'])

% Micro_coordinate : Coordinates of each mesh point
% Micro_element       : Mesh conectivity

global xi Macro_geo Micro_geo Time
% xi         : HDG Parameter

Macro_geo = struct();
Micro_geo = struct();

xi = 1;   % xi > 0

% coordinate ->  [ x , y ]
% Gamma      ->  Boundary  [ Label , NodeIni, NodeFin]
% element    ->  Element [ Node1 , Node2 , Node3]
[Macro_geo.coordinate, MacroGamma, Macro_geo.element] = readGmsh([pwd,meshfolder,file1]);
[Micro_geo.coordinate, MicroGamma, Micro_geo.element] = readGmsh([pwd,meshfolder,file2]);

Micro_geo.coordinate = Micro_geo.coordinate(:,1:2);
Macro_geo.coordinate = Macro_geo.coordinate(:,1:2);

% nElement -> number of element at each mesh
Macro_geo.nElement  = size(Macro_geo.element,1);
Micro_geo.nElement  = size(Micro_geo.element,1);

% nnodes -> number of nodes at each mesh
Macro_geo.nnodes = size(Macro_geo.coordinate,1);
Micro_geo.nnodes = size(Micro_geo.coordinate,1);

%*** ----------------------------------------------------------------
% TIME FEATURES
% dt : Time step discretization
% time_vect = t_0 : dt : T

Time.MicroTdepend = 1;
Time.dt           = 0.05;
Time.time_vec     = 0:Time.dt:1;
Time.tnSteps      = length(Time.time_vec);
Time.Inicial      = InicialSolution;

[Time.Inicial_cont,~] = PostProcess(Macro_geo,Time.Inicial,zeros(6*Macro_geo.nElement,1),1);

%%
%*********************************************************************
%*                                                                   *
%*                 Relation between Element and Edges                *
%*                                                                   *
%*********************************************************************

% (Paper: Carstensen - RaviarThomas elements)
% nodes2edge   : dim = nnodesxnnodes
% nodes2edge(i,j) = k where k is teh number of the edge between i,j.
%
% edge2element : dim = nedges x 4
% edge2element(i,:) -> [Inicial Final Element1 Element2 (0 Boundary)]
%
% interioredge : interioredge = edge2element only interior

% MACRO features
[~,Macro_geo.nodes2edge,~,Macro_geo.edge2element,Macro_geo.interioredge] = ...
    edge(Macro_geo.element,Macro_geo.coordinate);

% MICRO features
[~,Micro_geo.nodes2edge,~,Micro_geo.edge2element,Micro_geo.interioredge] = ...
    edge(Micro_geo.element,Micro_geo.coordinate);

%---------------------------------------------------------------------

% Sorting Boundaries

% NewmanEdges -DirichletEdges : dim = nedges_FronteraNewmann x 4
% NewmanEdges/DirichletEdges(i,1) -> [counter Initial Final Element]

% ############# this can be unusefull, is importante to review the way to
% ############# read/generate meshes and generate mesh structures.
% ############# This code can/should be optimizated.

[Micro_geo.NewmanEdges,Micro_geo.DirichletEdges] = ...
    sortBoundary(Micro_geo,MicroGamma,[indGamma_D2,indGamma_N2]);
[Macro_geo.NewmanEdges,Macro_geo.DirichletEdges] = ...
    sortBoundary(Macro_geo,MacroGamma,[indGamma_D1,indGamma_N1]);

Micro_geo.nEdgeInt    = size(Micro_geo.interioredge,1);
Micro_geo.nEdgeDir    = size(Micro_geo.DirichletEdges,1);
Micro_geo.nEdgeNew    = size(Micro_geo.NewmanEdges,1);
Micro_geo.nEdgeTotal  = Micro_geo.nEdgeInt+Micro_geo.nEdgeDir+Micro_geo.nEdgeNew;

Macro_geo.nEdgeInt    = size(Macro_geo.interioredge,1); %Number of int Edges
Macro_geo.nEdgeDir    = size(Macro_geo.DirichletEdges,1); %Number of dir Edges
Macro_geo.nEdgeNew    = size(Macro_geo.NewmanEdges,1); %Number of newmann Edges
Macro_geo.nEdgeTotal  = Macro_geo.nEdgeInt+Macro_geo.nEdgeDir+Macro_geo.nEdgeNew;

%*********************************************************************
% Identification of the micro Boundaries to impose the periodicity   *
% condition if it is necesary/request.                               *
%*********************************************************************

if strcmp(MicroConst,'Periodic')
    % Micro_geo.edgePar -> pairs of edges on opposite sides of the borders.
    
    Micro_geo = PreProcess_Identification(Micro_geo);
end

Micro_geo = Position_indicators(Micro_geo,0);
Macro_geo = Position_indicators(Macro_geo,1);

% Preprocess
% TimePreprocess = toc;
% fprintf('\n Preprocess Time %.2f seg \n',TimePreprocess)

%%
%*********************************************************************
%*                                                                   *
%*                            HDG SOLVER                             *
%*                                                                   *
%*********************************************************************
% addpath([cd,'\Inputs'])
addpath([cd,'\Solver'])

% Solution structure
Macro_Sol = struct();
field     = sprintf('time%i',0);
Macro_Sol.(field).Pres = Time.Inicial;
Macro_Sol.(field).Vel  = [];
Macro_Sol.(field).A_Efective  = 0;

% Auxiliar matrices _ independent of time
[Macro_Aux] = Macro_0time_mat;

if Time.MicroTdepend == 0
    [Micro_Sol] =  Micro_SolucionHDGP1P0([0,0],0);
    [A_Efective] = EfectivePermTensor([0,0],Micro_Sol,0);
end

% Analitical solution over the mesh
[Macro_solution,Macro_solution_x,Macro_solution_y,...
    ~, ~, ~,~, ~, ~] = solucionExacta;

% field     = sprintf('time%i',1);
for t_pos = 2:Time.tnSteps
    
    %     clc
    fprintf('\n TimeStep %i - %i\n',t_pos,Time.tnSteps)
    
    switch Adepend
        case 'Macro'
            
            % Previous solution
            Pre_ant = Macro_Sol.(field).Pres;
            Pre_Aefect = Macro_Sol.(field).A_Efective;
            field = sprintf('time%i',t_pos-1);
            
            [Macro_Sol.(field)] =  Macro_SolucionHDGP1P1...
                (Pre_ant,Macro_Aux,0,...
                Pre_Aefect,t_pos);
            
        case 'Micro'
            
            if Time.MicroTdepend == 1
                [Micro_Sol] =  Micro_SolucionHDGP1P0([0,0],tt);
                [A_Efective] = EfectivePermTensor([0,0],Micro_Sol,t_pos);
            end
            
            % Previous solution
            Pre_ant = Macro_Sol.(field).Pres;
            field = sprintf('time%i',t_pos-1);
            
            [Macro_Sol.(field)] =  Macro_SolucionHDGP1P1...
                (Pre_ant,Macro_Aux,A_Efective,0,t_pos);
            
            %             % Post Process Micro cell problem
            %             [Micro_Sol.Pres1_cont,...
            %                 Micro_Sol.Vel1_cont,...
            %                 Micro_Sol.MagVel1_cont] = PostProcess...
            %                 (Micro_geo,Micro_Sol.Pres1,Micro_Sol.Vel1,0);
            %
            %             [Micro_Sol.Pres2_cont,...
            %                 Micro_Sol.Vel2_cont,...
            %                 Micro_Sol.MagVel2_cont] = PostProcess...
            %                 (Micro_geo,Micro_Sol.Pres2,Micro_Sol.Vel2,0);
            %
            %             MacroP = [0,0];
            %             run Plot_MicroSol
    end
    
    
    %%
    %*********************************************************************
    %*                                                                   *
    %*                 PostProcess  +  Errror                            *
    %*                                                                   *
    %*********************************************************************
    % Post Process Macro problem
    [Macro_Sol.(field).Pres_cont,...
        Macro_Sol.(field).Vel_cont,...
        Macro_Sol.(field).MagVel_cont] = PostProcess...
        (Macro_geo,Macro_Sol.(field).Pres,Macro_Sol.(field).Vel,1);
    
    % Pressure error
    Error_p   = errorL2_Macro(t_pos,Macro_Sol.(field).Pres,...
        Macro_solution,Macro_geo);
    % Reference norm
    Error_ref = errorL2_Macro(t_pos,zeros(size(Macro_Sol.(field).Pres,1),1),...
        Macro_solution,Macro_geo);
    % Relative error
    L2MacroError_rel(t_pos-1) = sqrt(Error_p)/sqrt(Error_ref)
    
    Error_dp    = errorH1P1P1(t_pos,Macro_Sol.(field).Vel,...
        Macro_solution_x,Macro_solution_y,Macro_geo);
    Error_dref  = errorH1P1P1(t_pos,zeros(size(Macro_Sol.(field).Vel,1),1),...
        Macro_solution_x,Macro_solution_y,Macro_geo);
    % Relative error
    H1MacroError_rel(t_pos-1) = sqrt(Error_p+Error_dp)/sqrt(Error_ref+Error_dref)
    
    aaasav = sprintf('Solut_time%i.mat',t_pos);
    save(aaasav)
end
%% GRAPHICS
addpath([cd,'\Graphics'])
%
% run Plot_MicroSol
run Plot_MacroSol_Time