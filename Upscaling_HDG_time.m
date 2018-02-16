%*********************************************************************
%*                                                                   *
%*                      HDG - Upscaling                              *
%*                                                                   *
%*                 Manuela Bastidas  - 2017                          *
%*               Universiteit Hasselt - CMAT                         *
%*                                                                   *
%*********************************************************************
%
% Two-dimensional  HDG method for elliptic problems.
%
% HDG_UP solves the equation
%     - div(a*grad(u)) = f in Omega
%                    u = 0 on the Dirichlet boundary
%               d/dn u = 0   on the Neumann boundary
% with tensor a(x) on a geometry described by trian.
%
% See folder Inputs
%
%***------------------------------------
% Manuela Bastidas  - 2017.

%%
%*********************************************************************
%*                                                                   *
%*                      Simulation Parameters                        *
%*                                                                   *
%*********************************************************************

% Preprocess time start
% tic;

function [Macro_Sol, Errorinfo] = Upscaling_HDG_time(Adepend,Time,Macro_geo,Micro_geo) 


%%
%*********************************************************************
%*                                                                   *
%*                            HDG SOLVER                             *
%*                                                                   *
%*********************************************************************
% addpath([cd,'\Inputs'])
% addpath([cd,'\Solver'])

% Solution structure
Macro_Sol = struct();
field     = sprintf('time%i',0);
Macro_Sol.(field).Pres = Time.Inicial;
Macro_Sol.(field).Vel  = [];
Macro_Sol.(field).A_Efective  = 0;

% Auxiliar matrices _ independent of time
[Macro_Aux] = Macro_0time_mat(Macro_geo);

if Time.MicroTdepend == 0
    [Micro_Sol] =  Micro_SolucionHDGP1P0([0,0],0);
    [A_Efective] = EfectivePermTensor([0,0],Micro_Sol,0);
end

% Analitical solution over the mesh
[Macro_solution,Macro_solution_x,Macro_solution_y,...
    ~, ~, ~,~, ~, ~] = solucionExacta;

Error_p = zeros(Time.tnSteps,1);
Error_ref = zeros(Time.tnSteps,1);
Error_dp = zeros(Time.tnSteps,1);
Error_dref = zeros(Time.tnSteps,1);

% field     = sprintf('time%i',1);
for t_pos = 2:Time.tnSteps
    
    switch Adepend
        case 'Macro'
            
            % Previous solution
            Pre_ant = Macro_Sol.(field).Pres;
            Pre_Aefect = Macro_Sol.(field).A_Efective;
            field = sprintf('time%i',t_pos-1);
            
            [Macro_Sol.(field)] =  Macro_SolucionHDGP1P1...
                (Pre_ant,Macro_Aux,0,...
                Pre_Aefect,t_pos,Macro_geo,Adepend,Time,Micro_geo);
            
        case 'Micro'
            
            if Time.MicroTdepend == 1
                [Micro_Sol] =  Micro_SolucionHDGP1P0([0,0],tt);
                [A_Efective] = EfectivePermTensor([0,0],Micro_Sol,t_pos);
            end
            
            % Previous solution
            Pre_ant = Macro_Sol.(field).Pres;
            field = sprintf('time%i',t_pos-1);
            
            [Macro_Sol.(field)] =  Macro_SolucionHDGP1P1...
                (Pre_ant,Macro_Aux,A_Efective,0,t_pos,Macro_geo,Adepend,Time,Micro_geo);
            
            %             % Post Process Micro cell problem
            %             [Micro_Sol.Pres1_cont,...
            %                 Micro_Sol.Vel1_cont,...
            %                 Micro_Sol.MagVel1_cont] = PostProcess...
            %                 (Micro_geo,Micro_Sol.Pres1,Micro_Sol.Vel1,0);
            %
            %             [Micro_Sol.Pres2_cont,...
            %                 Micro_Sol.Vel2_cont,...
            %                 Micro_Sol.MagVel2_cont] = PostProcess...
            %                 (Micro_geo,Micro_Sol.Pres2,Micro_Sol.Vel2,0);
            %
            %             MacroP = [0,0];
            %             run Plot_MicroSol
    end
    
    
    %%
    %*********************************************************************
    %*                                                                   *
    %*                 PostProcess  +  Errror                            *
    %*                                                                   *
    %*********************************************************************
    % Post Process Macro problem
    [Macro_Sol.(field).Pres_cont,...
        Macro_Sol.(field).Vel_cont,...
        Macro_Sol.(field).MagVel_cont] = PostProcess...
        (Macro_geo,Macro_Sol.(field).Pres,Macro_Sol.(field).Vel,1);
    
    % Pressure error
    Error_p(t_pos)   = errorL2_Macro...
        (t_pos,Macro_Sol.(field).Pres,...
        Macro_solution,Macro_geo);
    % Reference norm
    Error_ref(t_pos) = errorL2_Macro...
        (t_pos,zeros(size(Macro_Sol.(field).Pres,1),1),...
        Macro_solution,Macro_geo);
    % Relative error
%     L2MacroError_rel(t_pos-1) = sqrt(Error_p(t_pos-1))/sqrt(Error_ref(t_pos-1));
    
    Error_dp(t_pos)    = errorH1P1P1(t_pos,Macro_Sol.(field).Vel,...
        Macro_solution_x,Macro_solution_y,Macro_geo);
    Error_dref(t_pos)  = errorH1P1P1(t_pos,zeros(size(Macro_Sol.(field).Vel,1),1),...
        Macro_solution_x,Macro_solution_y,Macro_geo);
    % Relative error
%     H1MacroError_rel(t_pos-1) = sqrt(Error_p(t_pos-1)+Error_dp(t_pos-1))/...
%         sqrt(Error_ref(t_pos-1)+Error_dref(t_pos-1));
    
%     aaasav = sprintf('Solut_time%i.mat',t_pos);
%     save(aaasav)
end

Errorinfo = {Error_p,Error_ref,Error_dp,Error_dref};
%% GRAPHICS
% addpath([cd,'\Graphics'])
% %
% % run Plot_MicroSol
% run Plot_MacroSol_Time