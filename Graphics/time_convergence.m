%% Convergence plot - Time to stationary 

pendiente = diff(L2MacroError_rel_inf)./diff(Time.time_vec(2:end))
figure 

plot(Time.time_vec(2:end),L2MacroError_rel_inf,'.-','linewidth',1)
hold on
% plot(time20,L2MacroError_rel_20tstep,'.-','linewidth',1)
% plot(Time.time_vec(2:end),H1MacroError_rel,'-','linewidth',1.5)
legend('10 steps error','20 steps error')
xlabel('Time step','FontSize',12,'FontWeight','bold')
ylabel('Error','FontSize',12,'FontWeight','bold')
% xtickformat('%.1f')
% ytickformat('%.1f')
grid on

%%
global Time
for t_pos = 2:Time.tnSteps
    field = sprintf('time%i',t_pos-1);
      % Pressure error
    Error_pL2(t_pos-1)   = errorL2_Macro(t_pos,Macro_Sol.(field).Pres,...
        Macro_solution,Macro_geo);
    % Reference norm
    Error_refL2(t_pos-1) = errorL2_Macro(t_pos,zeros(size(Macro_Sol.(field).Pres,1),1),...
        Macro_solution,Macro_geo);
    % Relative error
%     L2MacroError_rel_inf(t_pos-1) = Error_p_inf/Error_ref_inf
end
ERROR = sqrt(Time.dt*sum(Error_pL2./Error_refL2))
