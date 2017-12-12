%*********************************************************************
%*                                                                   *
%*     Post process of the HDG solution Micro - Macro Scale          *
%*                                                                   *
%*********************************************************************
%
% This code transforms the discontinuous solution into a continuous
% solution. (The result is usefull to visualize the results)
%
%***------------------------------------
%***Input:    Geometry, Scalar Discontinous solution, Vectorial
%             Discontinous solution.
%
%***------------------------------------
%***Output:   Scalar Solution (Continuous)
%             Vectorial Solution  (Continuous)
%
%***------------------------------------
% Manuela Bastidas  - 2017.

function [Pressure_cont,Vel_cont,MagVel_cont] = PostProcess(geo,Pressure,Velocity,grad)

nElement = geo.nElement;
nnodes   = geo.nnodes;
element  = geo.element;

%%
%*********************************************************************
%*                                                                   *
%*          Magnitud of the velocity / vectorial solution            *
%*                                                                   *
%*********************************************************************

MagVel_aux = zeros(nnodes,1);
if grad ==0
    for j=1:nElement
        posdux = 2*j-1;
        posduy = 2*j;
        MagVel_aux(j,1)  = sqrt(Velocity(posdux).^2+...
            Velocity(posduy).^2);
    end
else
    for j=1:nElement
        posdux = 6*j-5:6*j-3;
        posduy = 6*j-2:6*j;
        MagVel_aux(3*j-2:3*j,1)  = sqrt(Velocity(posdux).^2+...
            Velocity(posduy).^2);
    end
end

%%
%*********************************************************************
%*                                                                   *
%*              Averanging the solution of each node                 *
%*                                                                   *
%*********************************************************************

Pressure_cont = zeros(nnodes,1);
MagVel_cont   = zeros(nnodes,1);
Vel_cont       = zeros(nnodes,2);

contador = zeros(nnodes,1);
for j=1:nElement
    pos1 = 3*(j-1)+1:3*j;
    
    Pressure_cont(element(j,:))= Pressure_cont(element(j,:)) + Pressure(pos1);
    if grad == 0
        posdux = 2*j-1;
        posduy = 2*j;
        MagVel_cont(element(j,:))  = MagVel_cont(element(j,:)) + MagVel_aux(j);
        Vel_cont(element(j,:),1)   = Vel_cont(element(j,:),1) + Velocity(posdux);
        Vel_cont(element(j,:),2)   = Vel_cont(element(j,:),2) + Velocity(posduy);
    else
        posdux = 6*j-5:6*j-3;
        posduy = 6*j-2:6*j;
        MagVel_cont(element(j,:))  = MagVel_cont(element(j,:)) + MagVel_aux(pos1);
        Vel_cont(element(j,:),1)   = Vel_cont(element(j,:),1) + Velocity(posdux);
        Vel_cont(element(j,:),2)   = Vel_cont(element(j,:),2) + Velocity(posduy);
    end
    % Counter of the number of element for each node.
    contador(element(j,:))     = contador(element(j,:)) + ones(3,1);
end
Pressure_cont = Pressure_cont./contador;
MagVel_cont   = MagVel_cont./contador;
Vel_cont(:,1) =  Vel_cont(:,1)./contador;
Vel_cont(:,2) =  Vel_cont(:,2)./contador;
