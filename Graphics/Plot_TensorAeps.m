%*********************************************************************
%*                                                                   *
%*              Plot the original Difussion Tensor                   *
%*                                                                   *
%*********************************************************************
%
% This code is the auxiliary code to graph the behavior of the initial 
% oscillatory diffusion tensor.
%
% See function: TensorA_eps
%
%***------------------------------------
% Manuela Bastidas  - 2017.


%% Tensor evaluation 

% Macro Point
punto = [0,0];
% Mesh - Grid of [0,1]^2 
[X,Y] = meshgrid(0:0.005:1,0:0.005:1);

% Tensor evaluation
[TensorA,~,RHS] = TensorA_eps(punto,[X(:),Y(:)]);

% Components of the tensor
Z_11 = reshape(TensorA(1,1,:),size(X,1),size(Y,2));
Z_12 = reshape(TensorA(1,2,:),size(X,1),size(Y,2));
Z_21 = reshape(TensorA(2,1,:),size(X,1),size(Y,2));
Z_22 = reshape(TensorA(2,2,:),size(X,1),size(Y,2));


%% Plot

figure('name','Difussion Tensor')
surf(X,Y,Z_11,'Edgecolor','none')
xlabel('x1'); ylabel('x2')
view(0,90)
colorbar
colormap jet
title('Tensor A')


