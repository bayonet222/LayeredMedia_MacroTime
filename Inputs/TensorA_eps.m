%*********************************************************************
%*                                                                   *
%*             Difussion/Permeability Tensor (Highly oscilations)    *
%*                                                                   *
%*********************************************************************
%
% This function should be handle by the user to define the diffusion
% tensor that describes the oscilations in a micro scale.
%
% Remark: This function evaluate the tensor in several micro points.
%
%***------------------------------------
%***Input:    X: Macro Coordinate to eval the tensor
%             Y: Micro Coordinates to eval the tensor (Vector)
%
%***------------------------------------
%***Output:   TensorA     -> TensorA(:,:,i) = Tensor evaluated at x,y(i)
%             TensorA_inv -> Inverse Tensor.
%             RHS_cell    -> Right hand part of the cell problems [P1(i)
%             P2(i)]
%
%***------------------------------------
% Manuela Bastidas  - 2017.

function [TensorA,TensorA_inv,RHS_cell] = TensorA_eps(X,Y,tt)

% Each component of the tensor (Input of the original problem

% MICRO TEST 1
a_11 = @(x,y,tt) ( 1./(2-sin(2*pi*y(2))) ).*(exp(-tt)./(1+x(1).^2));
a_12 = @(x,y,tt) 0*y(1);
a_22 = @(x,y,tt) ( 1./(2-sin(2*pi*y(2))) ).*(exp(-tt)./(1+x(1).^2));
a_21 = @(x,y,tt) 0*y(1);

% Derivative respect y_j of A_ik
dy1_a11 = @(x,y,tt)0*y(1);
dy1_a12 = @(x,y,tt)0*y(1);

dy2_a22 = @(x,y,tt)( (2*pi.*cos(2*pi*y(2)))./((2-sin(2*pi*y(2))).^2)).*...
    (exp(-tt)./(1+x(1).^2));
% dy2_a22 = @(x,y)0*y(1);
dy2_a21 = @(x,y,tt)0*y(1);

%%  Inicialization

nPoints     = size(Y,1);
TensorA     = zeros(2,2,nPoints);
TensorA_inv = zeros(2,2,nPoints);
RHS_cell    = zeros(nPoints,2);

%% Evaluation

for i=1:nPoints
    
    %Point by Point
    y = Y(i,:);
    
    TensorA(:,:,i) = [a_11(X,y,tt),a_12(X,y,tt);
        a_21(X,y,tt),a_22(X,y,tt)];
    
    TensorA_inv(:,:,i) = (1./det(TensorA(:,:,i))).*...
        [a_22(X,y,tt),-a_12(X,y,tt);
        -a_21(X,y,tt),a_11(X,y,tt)];
    
    RHS_cell(i,:) = [dy1_a11(X,y,tt)+dy2_a21(X,y,tt),...
        dy1_a12(X,y,tt)+dy2_a22(X,y,tt)];
    
end

