%*********************************************************************
%*                                                                   *
%*              Matrix A of the macro local solver                   *
%*                                                                   *
%*********************************************************************
%
% This functions compute the matrix A (of HDG the local solvers)  for a
% macro problem. This is requiered because this matrix is the unique one
% that depends of the difussion tensor.
%
%               A(i,j) = int_K (phi_i) inv(Aperm*) (phi_j)
%
% Remark: We use P_0 aproximation for the vectorial solution of the
% problem, for this reason we will use only one point to integrate over the
% triangle.
%
%***------------------------------------
%***Inputs: coordinates of the triangle to compute the local solver
%
%***------------------------------------
% Manuela Bastidas - 2017.

function [matA_localsolver,A_Efective] = Macro_matA_localsolver2(Wx_macro,Wy_macro,...
    test,ind,MacroP,tt,Micro_geo)


% X_macro = X_macro(:);
% Y_macro = Y_macro(:);
nn = size(test,2);

[Micro_Sol] = Micro_SolucionHDGP1P0(MacroP,tt,Micro_geo);
   
[A_Efective] = EfectivePermTensor(MacroP,Micro_Sol,tt,Micro_geo);
inv_AEfective_mean = inv(A_Efective);

% A_Efective = repmat(A_Efective_mean,1,1,nn);
inv_AEfective = repmat(inv_AEfective_mean,1,1,nn);

%% MACRO - Numerical integration

if ind == 0
    matA_localsolver(1,1) = Wx_macro'*reshape(inv_AEfective(1,1,:),size(Wx_macro,1),size(Wy_macro,1))*Wy_macro;
    matA_localsolver(1,2) = Wx_macro'*reshape(inv_AEfective(2,1,:),size(Wx_macro,1),size(Wy_macro,1))*Wy_macro;
    matA_localsolver(2,1) = Wx_macro'*reshape(inv_AEfective(1,2,:),size(Wx_macro,1),size(Wy_macro,1))*Wy_macro;
    matA_localsolver(2,2) = Wx_macro'*reshape(inv_AEfective(2,2,:),size(Wx_macro,1),size(Wy_macro,1))*Wy_macro;
else
    
    for i=1:3
        for j=1:3
            multi = test(i,:).*test(j,:);
            
            uno = reshape(inv_AEfective(1,1,:),size(Wx_macro,1),size(Wy_macro,1)).*...
                reshape(multi,size(Wx_macro,1),size(Wy_macro,1));
            matA_localsolver1(i,j) = Wx_macro'*uno*Wy_macro;
            matA_localsolver1(j,i) = matA_localsolver1(i,j);
            
            uno = reshape(inv_AEfective(2,1,:),size(Wx_macro,1),size(Wy_macro,1)).*...
                reshape(multi,size(Wx_macro,1),size(Wy_macro,1));
            matA_localsolver2(i,j) = Wx_macro'*uno*Wy_macro;
            matA_localsolver2(j,i) = matA_localsolver2(i,j);
            
            uno = reshape(inv_AEfective(1,2,:),size(Wx_macro,1),size(Wy_macro,1)).*...
                reshape(multi,size(Wx_macro,1),size(Wy_macro,1));
            matA_localsolver3(i,j) = Wx_macro'*uno*Wy_macro;
            matA_localsolver3(j,i) = matA_localsolver3(i,j);
            
            uno = reshape(inv_AEfective(2,2,:),size(Wx_macro,1),size(Wy_macro,1)).*...
                reshape(multi,size(Wx_macro,1),size(Wy_macro,1));
            matA_localsolver4(i,j) = Wx_macro'*uno*Wy_macro;
            matA_localsolver4(j,i) = matA_localsolver4(i,j);
        end
    end
    matA_localsolver = [ matA_localsolver1 matA_localsolver2; matA_localsolver3 matA_localsolver4];
end


end