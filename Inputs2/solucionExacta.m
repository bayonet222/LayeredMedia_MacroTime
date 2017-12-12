%*********************************************************************
%*                                                                   *
%*                     Analytical solution
%*                                                                   *
%*********************************************************************
%
% This function defines te analytical solution of the pressure and the
% velocity field. This is useful to compute the error.
%
%***------------------------------------  
% Manuela Bastidas - 2017. 

function [Macro_sol,Macro_solx,Macro_soly,...
    Micro_sol1, Micro_solx1, Micro_soly1,...
    Micro_sol2, Micro_solx2, Micro_soly2] = solucionExacta

% sol =  @(x,y) x.*(x-1)+y.*(y-1);
% solx = @(x,y)2*x-1;
% soly = @(x,y) 2*y-1;

Macro_sol =  @(x1,x2,t) sin(pi*x1).*sin(pi*x2)+0.*t;
Macro_solx = @(x1,x2,t) pi*cos(pi*x1).*sin(pi*x2)+0.*t;
Macro_soly = @(x1,x2,t) pi*sin(pi*x1).*cos(pi*x2)+0.*t;

Micro_sol1 =  @(x1,x2,y1,y2) 0*x1+0*x2+0*y1+0*y2;
Micro_solx1 = @(x1,x2,y1,y2) 0*x1+0*x2+0*y1+0*y2;
Micro_soly1 = @(x1,x2,y1,y2) 0*x1+0*x2+0*y1+0*y2;

Micro_sol2 =  @(x1,x2,y1,y2) cos(2*pi*y2)./(4*pi)+0*x1+0*x2+0*y1;
Micro_solx2 = @(x1,x2,y1,y2) 0*x1+0*x2+0*y1+0*y2;
Micro_soly2 = @(x1,x2,y1,y2) -(sin(2*pi*y2))./2 +0*x1+0*x2+0*y1;



