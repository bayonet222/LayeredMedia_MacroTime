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

function [matA_localsolver,A_Efective_mean] = Macro_matA_localsolver(X_macro,Y_macro,Wx_macro,Wy_macro,...
    test,ind)

% global Micro_geo

% areaK = Macro_geo.areaK{e,1};
% N = 1 --> P_0
% Points and weight for the numerical integration.
% N=2;
% [X_macro,Y_macro,Wx_macro,Wy_macro] = triquad(N,coord');
X_macro = X_macro(:);
Y_macro = Y_macro(:);

A_Efective    = zeros(2,2,size(X_macro,1));
inv_AEfective = zeros(2,2,size(X_macro,1));

% Microerror = zeros(2,size(X_macro,1));
% Loop over the Macro integration points
for i = 1:size(X_macro,1)
    
    MacroP = [X_macro(i),Y_macro(i)];
    
    %% Micro cell problem for each Macro point
     fprintf('- %i.%i',i,size(X_macro,1))
    [Micro_Sol] = Micro_SolucionHDGP1P0(MacroP);
    
    %% A_efective computation
    % Integral over MICRO POROUS SPACE = sum(int_microelement .)
    [A_aux] = EfectivePermTensor(MacroP,Micro_Sol);

    A_Efective(:,:,i)    = A_aux;
    inv_AEfective(:,:,i) = inv(A_aux);
    
%     [~,~,~,Micro_solution1,~,~,Micro_solution2,~,~] = solucionExacta;
%     
%     % Post Process Micro cell problem
%     [Micro_Sol.Pres1_cont,...
%         Micro_Sol.Vel1_cont,...
%         Micro_Sol.MagVel1_cont] = PostProcess...
%         (Micro_geo,Micro_Sol.Pres1,Micro_Sol.Vel1,0);
%     
%     [Micro_Sol.Pres2_cont,...
%         Micro_Sol.Vel2_cont,...
%         Micro_Sol.MagVel2_cont] = PostProcess...
%         (Micro_geo,Micro_Sol.Pres2,Micro_Sol.Vel2,0);
    
%     Micro_Sol = movesolution(MacroP,Micro_geo,Micro_solution1,Micro_solution2,Micro_Sol);
    
    % ***----------------------------------------
    % Micro ERROR
    % Pressure error
%     Error_p_micro1   = errorL2(MacroP,Micro_Sol.Pres1,Micro_solution1,Micro_geo);
%     Error_p_micro2   = errorL2(MacroP,Micro_Sol.Pres2,Micro_solution2,Micro_geo);
%     % Reference norm
%     Error_ref_micro1 = errorL2(MacroP,zeros(size(Micro_Sol.Pres1,1),1),Micro_solution1,Micro_geo);
%     Error_ref_micro2 = errorL2(MacroP,zeros(size(Micro_Sol.Pres2,1),1),Micro_solution2,Micro_geo);
%     
%     % Relative error
%     Microerror(1,i) = Error_p_micro1/(Error_ref_micro1+eps)
%     Microerror(2,i) = Error_p_micro2/(Error_ref_micro2+eps)
    
%     addpath([cd,'\Graphics'])
%     %
%     run Plot_MicroSol
end
A_Efective_mean = sum(A_Efective,3)./size(X_macro,1);
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