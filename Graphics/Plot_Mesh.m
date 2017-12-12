%*********************************************************************
%*                                                                   *
%*                      HDG - Upscaling                              *
%*                    MESHES VISUALIZATION                           *
%*                                                                   *
%*********************************************************************
%
% This code is usefull to visualizate the main features of the micro 
% and macro meshes.
%
%***------------------------------------
% Manuela Bastidas  - 2017.


%% Macro mesh

global Macro_geo Micro_geo

figure('name', 'Macro - Scale')
trimesh(Macro_geo.element,...
    Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
    zeros(size(Macro_geo.coordinate,1),1))
grid off
view(0,90)

%% Micro mesh 

figure('name', 'Micro - Scale')
trimesh(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    zeros(size(Micro_geo.coordinate,1),1))
grid off
view(0,90)
