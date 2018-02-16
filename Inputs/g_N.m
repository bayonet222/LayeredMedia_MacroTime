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

val = 0;
 