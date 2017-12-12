clc
close all

global Time Macro_geo

v = VideoWriter('Solution.avi');
v.FrameRate = 5;
vid.FramesPerTrigger = 5; %as frame rate
open(v);

hFig = figure(1);

axis tight manual
set(gca,'nextplot','replacechildren');
set(gcf,'color','w');

for tt = 2:Time.tnSteps
    %     pause(0.5)
    clf
    
    field = sprintf('time%i',tt-1);
    trisurf(Macro_geo.element,...
        Macro_geo.coordinate(:,1),Macro_geo.coordinate(:,2),...
        Macro_Sol.(field).Pres_cont,...
        'facecolor','interp','Edgecolor','interp');
    grid off
    
    colormap parula
    c1 = colorbar;
    ylabel(c1,'Pressure','fontsize',12)
    view(0,90)
    aa = sprintf('Macro Solution - Time %.1f',Time.time_vec(tt));
    title(aa,'fontsize',14)
    
    frame = getframe(hFig);
    writeVideo(v,frame);
end

close(v);