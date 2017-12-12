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
% clear
% close all
% clc


% global  MicroConst
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

% Micro_coordinate : Coordinates of each mesh point
% Micro_element       : Mesh conectivity

% global xi Macro_geo Micro_geo
% xi         : HDG Parameter
function [Micro_geo,Macro_geo] = preUpscaling_HDG(file1,file2)

meshfolder  = '/Mesh/UnitSquare/';
nombres     = get_list_files([cd,meshfolder],'');
% nombres2 = nombres;
% Adepend    = 'Macro';
MicroConst = 'Periodic';
indGamma_N2  = NaN;
indGamma_D2  = 7;
indGamma_N1  = NaN;
    indGamma_D1  = 7;
    
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

% %*** ----------------------------------------------------------------
% % TIME FEATURES
% % dt : Time step discretization
% % time_vect = t_0 : dt : T
% 
% Time.MicroTdepend = 1;
% Time.dt           = 0.05;
% Time.time_vec     = 0:Time.dt:1;
% Time.tnSteps      = length(Time.time_vec);
% Time.Inicial      = InicialSolution;
% 
% [Time.Inicial_cont,~] = PostProcess(Macro_geo,Time.Inicial,zeros(6*Macro_geo.nElement,1),1);

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

