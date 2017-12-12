%*********************************************************************
%*                       Inicial Solution u_0                        *
%*********************************************************************
%
%This functions returns the evaluation of the initial function u_0 at each
%point of the mesh (x element), this means that the initial solution is a
%discontinuous function.
%
%***------------------------------------  
% Manuela Bastidas  - 2017. 

function [sol] = InicialSolution(Macro_geo)

% global Macro_geo 
nElement   = Macro_geo.nElement;
element    = Macro_geo.element;
coordinate = Macro_geo.coordinate;

%% Initial function
% User definied 

solt0 = @(x1,x2) sin(pi*x1).*sin(pi*x2);
%%
sol = zeros(3*nElement,1);

for j=1:nElement
    pos = 3*(j-1)+1:3*j;
    for i=1:3
        xx = coordinate(element(j,i),1);
        yy = coordinate(element(j,i),2);
        sol(pos(i),1) = solt0(xx,yy);
    end
end    

