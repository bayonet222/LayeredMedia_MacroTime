%*********************************************************************
%*                                                                   *
%*                      Newmann boundary condition                   *
%*                                                                   *
%*********************************************************************
%
% This function receives the coordinates x = [x; y] of one point 
% in the domain and the normal vector. 
% It evaluates the inner product between the newman boundary function 
% and the normal. (dot(gn,normal))
%
%***------------------------------------
%***Input:  x: Coordinates to eval.
%           n: normal vector.
%***------------------------------------  
% Manuela Bastidas - 2017. 

function val = g_N(x,n)
 xx = x(1);
 yy = x(2);
 %% Problema Homogeneo
%  val = [yy*(yy-1)*((xx-1)+xx)  xx*(xx-1)*((yy-1)+yy)]*n';

 %% Problema  No Homogeneo
%  val = [(2*xx-1)  (2*yy-1)]*n';
  val1 = ([pi*cos(pi*xx).*sin(pi*yy); pi*sin(pi*xx).*cos(pi*yy)]);
  val = dot(val1,n);
% [pi*cos(pi*xx).*sin(pi*yy) pi*sin(pi*xx).*cos(pi*yy)];
%  val = [1 1]*n';
% val = 0;
 