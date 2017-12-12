%*********************************************************************
%*                                                                   *
%*                      HDG - Upscaling                              *
%*          Solution of the Macro Scale Upscaled problem             *
%*                                                                   *
%*********************************************************************
%
%***Main file: This code is the main to plot the results of the 
%              Macro scale problem.
%
% This code uses de results of the HDG_Upscaling.m and Pros-processing
%
%***------------------------------------
% Manuela Bastidas  - 2017.

global Macro_geo

%%
%*********************************************************************
%*                                                                   *
%*                         Scalar solution                           *
%*                                                                   *
%*********************************************************************
solZ = Macro_solution(Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2));

colormin = min(Macro_Sol.Pres_cont);
colormax = max(Macro_Sol.Pres_cont);

figure('name','Macro solution pressure')
h1 = subplot(1,2,1);
trisurf(Macro_geo.element,...
    Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
    Macro_Sol.Pres_cont,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap parula
caxis(h1,[colormin colormax])
c1 = colorbar;
ylabel(c1,'Pressure')
% ytickformat('%.1f')
% xtickformat('%.1f')
view(0,90)
title({'Macro problem';'HDG solution'})

h2=subplot(1,2,2);
trisurf(Macro_geo.element,...
    Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
    abs(Macro_Sol.Pres_cont-solZ),...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap parula
% caxis(h2,[colormin colormax])
c1 = colorbar;
ylabel(c1,'Pressure')
% ytickformat('%.1f')
% xtickformat('%.1f')
view(0,90)
title('Error')



%%
%*********************************************************************
%*                                                                   *
%*                         Vectorial solution                        *
%*                                                                   *
%*********************************************************************
solZx1 = Macro_solution_x(Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2));
solZy1 = Macro_solution_y(Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2));

SolMagxy_macro = abs(sqrt(solZx1.^2+solZy1.^2));

colormin = min(Macro_Sol.MagVel_cont);
colormax = max(Macro_Sol.MagVel_cont);

figure
h1 = subplot(1,2,1);
trisurf(Macro_geo.element,...
    Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
    Macro_Sol.MagVel_cont,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap parula
caxis(h1,[colormin colormax])
c1 = colorbar;
ylabel(c1,'$\| \nabla Pressure \|$','interpreter','latex')
% ylabel(c1,'Pressure Micro1','fontsize',16)
% ytickformat('%.1f')
% xtickformat('%.1f')
view(0,90)
title({'Macro problem';'HDG solution'})

h2=subplot(1,2,2);
trisurf(Macro_geo.element,...
    Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
    abs(Macro_Sol.MagVel_cont-SolMagxy_macro),...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap parula
% caxis(h2,[colormin colormax])
c1 = colorbar;
ylabel(c1,'$\| \nabla Pressure \|$','interpreter','latex')
% ylabel(c1,'Pressure Micro2','fontsize',16)
% ytickformat('%.1f')
% xtickformat('%.1f')
view(0,90)
title('Error')
