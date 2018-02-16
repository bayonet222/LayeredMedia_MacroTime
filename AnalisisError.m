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
% clear
clc

%% Set of meshes
%
% global Adepend MicroConst
% global xi Macro_geo Micro_geo

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

Solution = cell(4,1);
ErrorSol = cell(4,1);

for pqr = 6:10
    
    file1 = nombres{pqr};
    indGamma_N1  = NaN;
    indGamma_D1  = 7;
    
    fprintf('\n --------------------- \n')
    disp(file1)
    
    [Micro_geo,Macro_geo] = preUpscaling_HDG(file1,file2);
    
    % TIME DEPENDENCY
    Time = struct();
    Time.MicroTdepend = 1;
    Time.tnSteps      = ceil(1/Macro_geo.size);
    Time.dt           = 1/Time.tnSteps;
    Time.time_vec     = 0:Time.dt:1;
    Time.Inicial      = InicialSolution(Macro_geo);
    
    %% Error process
    [Macro_Sol, Errorinfo] = Upscaling_HDG_time(Adepend,...
        Time,Macro_geo,Micro_geo);
    
    Solution{pqr} = Macro_Sol;
    ErrorSol{pqr} = Errorinfo;
    
    aaa = sprintf('sol_%i.mat',pqr);
    save(aaa) 
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

