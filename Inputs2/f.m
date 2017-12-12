%*********************************************************************
%*                                                                   *
%*                      Flux - Source (f)                            *
%*                                                                   *
%*********************************************************************
%
% This function receives the coordinates x = [x; y] in the domain and
% evaluates the source associated with the equation to solve.
%
%***------------------------------------
%***Input:  x: Coordenadinates to eval.
%
%***------------------------------------  
% Manuela Bastidas - 2017. 

function volumeforce = f(xx,yy,tt)
% xx = x(1);
% yy = x(2);

% volumeforce = ((2+sqrt(3))/(2*sqrt(3)))*pi^2*sin(pi*xx).*sin(pi*yy);

term = xx.^2+1;
c1 = (1/sqrt(3)+1/2)*((pi.^2)/term);
c2 = (2*pi*xx)./(sqrt(3)*term.^2);
volumeforce = c1.*sin(pi*xx).*sin(pi*yy)+...
    c2.*sin(pi*yy).*cos(pi*xx);