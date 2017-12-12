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

global Macro_geo Time

%%
%*********************************************************************
%*                                                                   *
%*                         Scalar solution                           *
%*                                                                   *
%*********************************************************************
% solZ = Macro_solution(Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2));

field = sprintf('time%i',0);
colormin = min(Macro_Sol.(field).Pres);
colormax = max(Macro_Sol.(field).Pres);
for tt = 2:Time.tnSteps
% for tt = 2:3
    field = sprintf('time%i',tt-1);
    colormin = min(min(Macro_Sol.(field).Pres_cont),colormin);
    colormax = max(max(Macro_Sol.(field).Pres_cont),colormax);
end

% solZ = Macro_solution(Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2));

figure
for tt = 2:Time.tnSteps
%   for tt = 2:3  
    pause(0.1)
    clf
    
    field = sprintf('time%i',tt-1);
    trisurf(Macro_geo.element,...
        Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
        Macro_Sol.(field).Pres_cont,...
        'facecolor','interp','Edgecolor','interp');
    
    grid off
    colormap parula
    caxis([colormin colormax])
    c1 = colorbar;
    ylabel(c1,'Pressure')
    % ytickformat('%.1f')
    % xtickformat('%.1f')
    view(0,90)
    aa = sprintf('Macro problem - Time %.1f',Time.time_vec(tt));
    title({aa;'HDG solution'})
    
end


%%
%*********************************************************************
%*                                                                   *
%*                         Vectorial solution                        *
%*                                                                   *
%*********************************************************************
field = sprintf('time%i',0);
colormin = 0;
colormax = 0;
for tt = 2:Time.tnSteps
    field = sprintf('time%i',tt-1);
    colormin = min(min(Macro_Sol.(field).MagVel_cont),colormin);
    colormax = max(max(Macro_Sol.(field).MagVel_cont),colormax);
end

figure
for tt = 2:Time.tnSteps
    
    pause(0.1)
    clf
    
    field = sprintf('time%i',tt-1);
    trisurf(Macro_geo.element,...
        Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
        Macro_Sol.(field).MagVel_cont,...
        'facecolor','interp','Edgecolor','interp');
    
    grid off
    colormap parula
    caxis([colormin colormax])
    c1 = colorbar;
    ylabel(c1,'$\| \nabla Pressure \|$','interpreter','latex')
    % ytickformat('%.1f')
    % xtickformat('%.1f')
    view(0,90)
    aa = sprintf('Macro problem - Time %.1f',Time.time_vec(tt));
    title({aa;'HDG solution'})
    
end