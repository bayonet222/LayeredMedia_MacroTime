%*********************************************************************
%*                                                                   *
%*                      HDG - Upscaling                              *
%*              Solution of the Micro Scale problems                 *
%*                                                                   *
%*********************************************************************
%
%***Main file: This code is the main to plot the results of the 
%              Micro scale problems.
%
% This code uses de results of the HDG_Upscaling.m and Pros-processing
%
%***------------------------------------
% Manuela Bastidas  - 2017.

global Micro_geo 

%%
%*********************************************************************
%*                                                                   *
%*                         Vectorial solution                           *
%*                                                                   *
%*********************************************************************

figure
subplot(1,2,1)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    Micro_Sol.Vel1_cont(:,1),...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro1_x','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 1','fontsize',14)

subplot(1,2,2)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    Micro_Sol.Vel1_cont(:,2),...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro1_y','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 1','fontsize',14)

figure
subplot(1,2,1)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    Micro_Sol.Vel2_cont(:,1),...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro2_x','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 2','fontsize',14)

subplot(1,2,2)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    Micro_Sol.Vel2_cont(:,2),...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro2_y','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 2','fontsize',14)


solZx1 = Micro_solx1(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));
solZy1 = Micro_soly1(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));

solZx2 = Micro_solx2(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));
solZy2 = Micro_soly2(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));

figure
subplot(1,2,1)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    solZx1,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro1_x Exacta','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 1','fontsize',14)

subplot(1,2,2)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    solZy1,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro1_y Exacta','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 1','fontsize',14)

figure
subplot(1,2,1)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    solZx2,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro2_x Exacta','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 2','fontsize',14)

subplot(1,2,2)
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    solZy2,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap jet
c1 = colorbar;
ylabel(c1,'Pressure Micro2_y Exacta','fontsize',16)
ytickformat('%.1f')
xtickformat('%.1f')
view(0,90)
title('Micro Problem 2','fontsize',14)