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
%*                         Scalar solution                           *
%*                                                                   *
%*********************************************************************
solZ1 = Micro_solution1(repmat(MacroP(1),size(Micro_geo.coordinate(:,1),1),1),...
    repmat(MacroP(2),size(Micro_geo.coordinate(:,1),1),1),...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));
solZ2 = Micro_solution2(repmat(MacroP(1),size(Micro_geo.coordinate(:,1),1),1),...
    repmat(MacroP(2),size(Micro_geo.coordinate(:,1),1),1),...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));

colormin = min(min([Micro_Sol.Pres1_cont Micro_Sol.Pres2_cont]));
colormax = max(max([Micro_Sol.Pres1_cont Micro_Sol.Pres2_cont]));

figure('name','Micro solution pressure')
h1 = subplot(2,2,1);
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    Micro_Sol.Pres1_cont,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap parula
caxis(h1,[colormin colormax])
c1 = colorbar;
% ylabel(c1,'Pressure Micro1','fontsize',16)
% ytickformat('%.1f')
% xtickformat('%.1f')
view(0,90)
title({'Micro problem 1';'HDG solution'})

h2 = subplot(2,2,2);
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    Micro_Sol.Pres2_cont,...
    'facecolor','interp','Edgecolor','interp');
grid off
colormap parula
caxis(h2,[colormin colormax])
c1 = colorbar;
ylabel(c1,'Pressure')
% ytickformat('%.1f')
% xtickformat('%.1f')
view(0,90)
title({'Micro problem 2';'HDG solution'})

h3=subplot(2,2,3);
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    abs(Micro_Sol.Pres1_cont-solZ1),...
    'facecolor','interp','Edgecolor','interp');
%     abs(Micro_Sol.Pres1_cont-solZ1),...
grid off
colormap parula
% caxis(h3,[colormin colormax])
c1 = colorbar;
% ylabel(c1,'Pressure Micro2','fontsize',16)
% ytickformat('%.1f')
% xtickformat('%.1f')
view(0,90)
title('Error')

h4 =subplot(2,2,4);
trisurf(Micro_geo.element,...
    Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
    abs(Micro_Sol.Pres2_cont-solZ2),...
    'facecolor','interp','Edgecolor','interp');
%     abs(Micro_Sol.Pres2_cont-solZ2),...
grid off
colormap parula
% caxis(h4,[colormin colormax])
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

% solZx1 = Micro_solution1_x(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));
% solZy1 = Micro_solution1_y(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));
% 
% SolMagxy_micro1 = abs(sqrt(solZx1.^2+solZy1.^2));
% 
% solZx2 = Micro_solution2_x(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));
% solZy2 = Micro_solution2_y(Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2));
% 
% SolMagxy_micro2 = (sqrt(solZx2.^2+solZy2.^2));
% 
% colormin = min(min([Micro_Sol.MagVel1_cont Micro_Sol.MagVel2_cont]));
% colormax = max(max([Micro_Sol.MagVel1_cont Micro_Sol.MagVel2_cont]));
% 
% figure('name','Micro solution gradient')
% h1 = subplot(2,2,1);
% trisurf(Micro_geo.element,...
%     Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
%     Micro_Sol.MagVel1_cont,...
%     'facecolor','interp','Edgecolor','interp');
% grid off
% colormap parula
% caxis(h1,[colormin colormax])
% c1 = colorbar;
% % ylabel(c1,'Pressure Micro1','fontsize',16)
% % ytickformat('%.1f')
% % xtickformat('%.1f')
% view(0,90)
% title({'Micro problem 1';'HDG solution'})
% 
% h2 = subplot(2,2,2);
% trisurf(Micro_geo.element,...
%     Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
%     Micro_Sol.MagVel2_cont,...
%     'facecolor','interp','Edgecolor','interp');
% grid off
% colormap parula
% caxis(h2,[colormin colormax])
% c1 = colorbar;
% ylabel(c1,'$\| \nabla Pressure \|$','interpreter','latex')
% % ytickformat('%.1f')
% % xtickformat('%.1f')
% view(0,90)
% title({'Micro problem 2';'HDG solution'})
% 
% h3=subplot(2,2,3);
% trisurf(Micro_geo.element,...
%     Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
%     abs(Micro_Sol.MagVel1_cont-SolMagxy_micro1),...
%     'facecolor','interp','Edgecolor','interp');
% grid off
% colormap parula
% % caxis(h3,[colormin colormax])
% c1 = colorbar;
% % ylabel(c1,'Pressure Micro2','fontsize',16)
% % ytickformat('%.1f')
% % xtickformat('%.1f')
% view(0,90)
% title('Error')
% 
% h4 =subplot(2,2,4);
% trisurf(Micro_geo.element,...
%     Micro_geo.coordinate(:,1),Micro_geo.coordinate(:,2),...
%     abs(Micro_Sol.MagVel2_cont-SolMagxy_micro2),...
%     'facecolor','interp','Edgecolor','interp');
% grid off
% colormap parula
% % caxis(h4,[colormin colormax])
% c1 = colorbar;
% ylabel(c1,'$\| \nabla Pressure \|$','interpreter','latex')
% % ytickformat('%.1f')
% % xtickformat('%.1f')
% view(0,90)
% title('Error')