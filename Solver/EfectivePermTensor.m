%*********************************************************************
%*                                                                   *
%*              EFECTIVE Diffusion tensor A _ Upscaled               *
%*                                                                   *
%*********************************************************************
%
% This functions compute the efective difussion tensor if the original
% tensor does not dependent of the macro scale.
%
%***------------------------------------
%***Inputs: Solution of the micro cell problems. 
%
%***------------------------------------
% Manuela Bastidas - 2017.

function [A_Efective] = EfectivePermTensor(MacroP,Micro_Sol,tt)

global Micro_geo 

Vel1 = Micro_Sol.Vel1;
Vel2 = Micro_Sol.Vel2;

AreaSum    = 0;
A_Efective = zeros(2,2);

for elemento = 1:Micro_geo.nElement
    
    % Coord (x;y) triangle vertex
    coord =  Micro_geo.coordinate( Micro_geo.element(elemento,:),:)';
    % Vectorial soluton position
    pos   = 2*elemento-1:2*elemento;
    
    % Integration points and weights. Quadrature rule
    [X_int,Y_int,Wx,Wy] = triquad(4,coord');
    [TensorA,~,~]       = TensorA_eps(MacroP,[X_int(:),Y_int(:)],tt);
    
    % A_EFECTIVE = INT_P (e_j+w^j)*e_i 
    AreaSum = AreaSum + Micro_geo.areaK{elemento,1};
    
    int_Aux11 = (Wx'*reshape(TensorA(1,1,:),size(Wx,1),size(Wy,1))*Wy);
    int_Aux12 = (Wx'*reshape(TensorA(1,2,:),size(Wx,1),size(Wy,1))*Wy);
    Aux(1,1) = [int_Aux11 int_Aux12]*(Vel1(pos)+[1;0]);
    Aux(1,2) = [int_Aux11 int_Aux12]*(Vel2(pos)+[0;1]);
    
    int_Aux21 = (Wx'*reshape(TensorA(2,1,:),size(Wx,1),size(Wy,1))*Wy);
    int_Aux22 = (Wx'*reshape(TensorA(2,2,:),size(Wx,1),size(Wy,1))*Wy);
    Aux(2,1) = [int_Aux21 int_Aux22]*(Vel1(pos)+[1;0]);
    Aux(2,2) = [int_Aux21 int_Aux22]*(Vel2(pos)+[0;1]);
    
    A_Efective = A_Efective + Aux;
end
A_Efective =  A_Efective./AreaSum;