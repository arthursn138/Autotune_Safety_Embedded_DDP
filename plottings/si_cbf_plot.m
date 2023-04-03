function si_cbf_plot(X, Xcbf, X_ddp, U, Ucbf, U_ddp, T, x0, xf, circ, filter, w_back, w_fwd, obs_loc)

IC_leg = ['IC: (' num2str(x0(1)) ',' num2str(x0(2)) '), Coefficient for \alpha in QP = 0.05; \alpha Filter = 7'];
obs_leg = ['Center at (' num2str(obs_loc(1)) ', ' num2str(obs_loc(2)) '), radius = ' num2str(obs_loc(3)) ' and taget at (' num2str(xf(1)) ',' num2str(xf(2)) ')'];

figure(1)
% % % rectangle('Position',circ(1,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); hold on
% % % invisible = line(NaN, NaN, 'LineWidth', 2, 'LineStyle', '-', 'Color',[0.6350 0.0780 0.1840],'DisplayName', obs_leg); % Workaround to have obstacle's info displayed
% % % plot(x0(1),x0(2),'*','LineWidth',1.5,'DisplayName','Initial State'); grid on;
% % % plot(xf(1),xf(2),'*','LineWidth',1.5,'DisplayName', 'Target');

% % % % % % % % % % % % % % % % plot(x0(1),x0(2),'*','LineWidth',1.5,'HandleVisibility','off'); grid on;
% % % % % % % % % % % % % % % % plot(xf(1),xf(2),'*','LineWidth',1.5,'HandleVisibility','off');
% % % % % % % % % % % % % % % % plot(X(1,:),X(2,:),'LineWidth',1.5,'DisplayName', IC_leg); 
% % % % % % % % % % % % % % % % % % % plot(X_ddp(1,:),X_ddp(2,:),':','LineWidth',1.5);
% % % % % % % % % % % % % % % % title('X-Y path (top view)','Interpreter','latex');
% % % % % % % % % % % % % % % % ylabel('$y$','FontName','Times New Roman', 'Interpreter','latex');
% % % % % % % % % % % % % % % % xlabel('$x$','FontName','Times New Roman', 'Interpreter','latex');
% % % % % % % % % % % % % % % % if filter
% % % % % % % % % % % % % % % %     plot(Xcbf(1,:),Xcbf(2,:),'--','LineWidth',1.5);
% % % % % % % % % % % % % % % %     legend('Initial state', 'Desired final state', 'CBF-Aware DDP', 'Vanilla DPP',...
% % % % % % % % % % % % % % % %         'CBF Filtering', 'FontName','Times New Roman','Interpreter','latex','Location','best');
% % % % % % % % % % % % % % % % else
% % % % % % % % % % % % % % % %     % % %     legend('Initial state', 'Desired final state', 'CBF-Aware DDP Path',...
% % % % % % % % % % % % % % % %     % % %         'Vanilla DPP Path', 'FontName','Times New Roman','Interpreter','latex');
% % % % % % % % % % % % % % % %     legend
% % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % box on; axis square; axis equal;
% % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % figure(2)
% % % % % % % % % % % % % % % % % % % % plot(T(1:end-1),U(1,:),'b','LineWidth',1.5); hold on; grid on;
% % % % % % % % % % % % % % % % % % % % plot(T(1:end-1),U(2,:),'r','LineWidth',1.5);
% % % % % % % % % % % % % % % % plot(T(1:end-1),U(1,:),'LineWidth',1.5,'DisplayName', ['u1,' IC_leg]); hold on; grid on;
% % % % % % % % % % % % % % % % plot(T(1:end-1),U(2,:),'LineWidth',1.5,'DisplayName', ['u2,' IC_leg]);
% % % % % % % % % % % % % % % % % % % plot(T(1:end-1),U_ddp(1,:),':b','LineWidth',1.5);
% % % % % % % % % % % % % % % % % % % plot(T(1:end-1),U_ddp(2,:),':r','LineWidth',1.5);
% % % % % % % % % % % % % % % % title('Control vs time','FontName','Times New Roman','Interpreter','latex');
% % % % % % % % % % % % % % % % ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
% % % % % % % % % % % % % % % % xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
% % % % % % % % % % % % % % % % if filter
% % % % % % % % % % % % % % % %     plot(T(1:end-1),Ucbf(1,:),'--b','LineWidth',1.5);
% % % % % % % % % % % % % % % %     plot(T(1:end-1),Ucbf(2,:),'--r','LineWidth',1.5);
% % % % % % % % % % % % % % % %     legend('$u_{1}$ CBF-Aware DDP', '$u_{2}$ CBF-Aware DDP', '$u_{1}$ vanilla',...
% % % % % % % % % % % % % % % %         '$u_{2}$ vanilla', '$u_{1}$ CBF Filtering', '$u_{2}$ CBF Filtering',...
% % % % % % % % % % % % % % % %         'FontName', 'Times New Roman', 'Interpreter','latex');
% % % % % % % % % % % % % % % % else
% % % % % % % % % % % % % % % % % % % %     legend('$u_{1}$ CBF-Aware', '$u_{2}$ CBF-Aware', '$u_{1}$ Vanilla DDP',...
% % % % % % % % % % % % % % % % % % % %         '$u_{2}$ Vanilla DDP', 'FontName','Times New Roman','Interpreter','latex');
% % % % % % % % % % % % % % % % % % % %     legend('$u_{1}$ CBF-Aware', '$u_{2}$ CBF-Aware', 'FontName','Times New Roman','Interpreter','latex');
% % % % % % % % % % % % % % % %     legend
% % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % box on; axis square; axis equal;


plot(X(1,:),X(2,:),'LineWidth',1.5,'DisplayName', ['CBF-Aware DDP: ' IC_leg]);
% % % % plot(X_ddp(1,:),X_ddp(2,:),':','LineWidth',1.5);
title('X-Y path (top view) [Min Norm]','Interpreter','latex');
ylabel('$y$','FontName','Times New Roman', 'Interpreter','latex');
xlabel('$x$','FontName','Times New Roman', 'Interpreter','latex');
plot(Xcbf(1,:),Xcbf(2,:),'--','LineWidth',1.5,'DisplayName', ['CBF Filtering: ' IC_leg]);
legend
box on; axis square; axis equal;

figure(2)
plot(T(1:end-1),U(1,:),'LineWidth',1.5,'DisplayName', ['u_1 CBF-Aware DDP ' IC_leg]); hold on; grid on;
plot(T(1:end-1),U(2,:),'LineWidth',1.5,'DisplayName', ['u_2 CBF-Aware DDP: ' IC_leg]);
% % % plot(T(1:end-1),U(1,:),'LineWidth',1.5,'DisplayName', ['u1,' IC_leg]); hold on; grid on;
% % % plot(T(1:end-1),U(2,:),'LineWidth',1.5,'DisplayName', ['u2,' IC_leg]);
% % % % % plot(T(1:end-1),U_ddp(1,:),':','LineWidth',1.5);
% % % % % plot(T(1:end-1),U_ddp(2,:),':','LineWidth',1.5);
title('Control vs time [Min Norm]','FontName','Times New Roman','Interpreter','latex');
ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
plot(T(1:end-1),Ucbf(1,:),'--','LineWidth',1.5,'DisplayName', ['u_1 CBF Filtering: ' IC_leg]);
plot(T(1:end-1),Ucbf(2,:),'--','LineWidth',1.5,'DisplayName', ['u_2 CBF Filtering: ' IC_leg]);
legend('Location','eastoutside')

figure(3)
plot(T, w_back,'LineWidth', 1.5,'DisplayName', IC_leg); hold on; grid on
title('$\alpha$ Back [Minimum Norm]','FontName','Times New Roman','Interpreter','latex'); 
ylabel('$\alpha$ back','FontName','Times New Roman','Interpreter','latex'); 
xlabel('time (s)','FontName','Times New Roman','Interpreter','latex');
legend

figure(4)
plot(T(1:end-1), w_fwd,'LineWidth', 1.5,'DisplayName', IC_leg); hold on; grid on
title('$\alpha$ Forward  [Minimum Norm]','FontName','Times New Roman','Interpreter','latex'); 
ylabel('$\alpha$ fwd','FontName','Times New Roman','Interpreter','latex'); 
xlabel('time (s)','FontName','Times New Roman','Interpreter','latex');
legend

end

