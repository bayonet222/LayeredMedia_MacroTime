%*********************************************************************
%*                                                                   *
%*              Error Analysis - Micro Problem HDG               *
%*                                                                   *
%*                 Manuela Bastidas  - 2017                          *
%*               Universiteit Hasselt - CMAT                         *
%*                                                                   *
%*********************************************************************
%
% Two-dimensional heterogeneous multiscale HDG method
% for elliptic problems.
%
% HDG_UP solves the equation
%     - div(a*grad(u)) = f in Omega
%                    u = 0 on the Dirichlet boundary
%               d/dn u = 0   on the Neumann boundary
% with tensor a(x) on a MICRO geometry described by trian.
%
% See folder Inputs
%
%***------------------------------------
% Manuela Bastidas  - 2017.
close all
clear
clc

%% Set of meshes
%
global Adepend MicroConst Time

% addpath([cd,'\Mesh'])
% 
% meshfolder  = '\Mesh\UnitSquare\';

addpath([cd,'/Mesh'])
addpath([cd,'/PrePost_process'])
addpath([cd,'/Inputs'])
addpath([cd,'/Solver'])

meshfolder  = '/Mesh/UnitSquare/';
nombres     = get_list_files([cd,meshfolder],'');
nombres2 = nombres;

%% Tabla -> Table nElement | mesh size | L2 error for each mesh
Tabla = zeros(size(nombres,1),4);

file2 = nombres{9};
Adepend    = 'Macro';
MicroConst = 'Periodic';
indGamma_N2  = NaN;
indGamma_D2  = 7;

Time.MicroTdepend = 1;

for pqr = 6:10
    
    file1 = nombres{pqr};
    indGamma_N1  = NaN;
    indGamma_D1  = 7;
    
    fprintf('\n --------------------- \n')
    disp(file1)
    
    % MESH PROCESS
    run preUpscaling_HDG
    
    % TIME DEPENDENCY
    Time.tnSteps      = ceil(1/Macro_geo.size);
    Time.dt           = 1/Time.tnSteps;
    Time.time_vec     = 0:Time.dt:1;
    Time.Inicial      = InicialSolution;
    
    %% Error process
    run Upscaling_HDG_time
    
    Tabla(pqr,1) = Macro_geo.nElement;
    Tabla(pqr,2) = Macro_geo.size;
    Tabla(pqr,3) = Time.dt;
    
    TablaError{pqr,1} = Error_p;
    TablaError{pqr,2} = Error_ref;
    Tabla(pqr,4) = sqrt(Time.dt*sum(Error_p./Error_ref));
    
    aaasav = sprintf('Solution_tot%i.mat',pqr);
    save(aaasav)
end


%% Convergence Report

% Tabla2 = sortrows(Tabla,1);
% Tabla2 = Tabla2(Tabla2(:,3)>0,:);
% 
% figure('name','L2 Convergence')
% hold on
% 
% axx = log(1./(Tabla2(:,2)+Tabla2(:,3)));
% plot(axx,...
%     log(Tabla2(:,4)),'s--','linewidth',1.5)
% % REFERENCE
% plot(axx, -1*axx-.55*min(abs(axx)),'linewidth',1.5)
% % plot(log(1./Tabla2(:,2)), -2*log(1./Tabla2(:,2))-.1*min(abs(log(Tabla2(:,3)))),'linewidth',1.5)
% 
% xlabel('Log(1/h)','FontSize',12,'FontWeight','bold')
% ylabel('log(Error)','FontSize',12,'FontWeight','bold')
% xtickformat('%.1f')
% ytickformat('%.1f')
% grid on
% xlim([min(axx)-0.05,max(axx)+0.05])

