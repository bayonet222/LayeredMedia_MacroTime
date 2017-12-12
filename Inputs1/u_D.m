%*********************************************************************
%*                                                                   *
%*                      Dirichlet boundary condition                 *
%*                                                                   *
%*********************************************************************
%
% This function receives the coordinates x = [x; y] of one point 
% in the domain and evaluates dirichlet function in that point.
%
%***------------------------------------
%***Input:  x: Coordinates to eval.
%
%***------------------------------------  
% Manuela Bastidas - 2017. 

function dir = u_D(x)
xx = x(1);
yy = x(2);

% dir =  xx.*(xx-1)+yy.*(yy-1);
% dir =   -(sin(pi*yy).^2)./4;
dir = 0;

