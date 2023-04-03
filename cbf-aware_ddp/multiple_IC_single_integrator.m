%% Control-Limited Differential Dynamic Programming
% Arthur Nascimento, Hassan Almubarak
% nascimento@gatech.edu, halmubarak@gatech.edu
% ACDS Lab @ Georgia Tech
% Last Update June 23 2022

% This is an example for 2D single integrator.
% This example uses optimal decay to allow for a better choice of alpha.

% CHECK BEFORE RUNNING:
% 0) Rename table to be saved (line 204)
% 1) Obstacle's location and size
% 2) size_of_set and permited range (factor that multiplies rand)
% 3) alpha, Q 
% 4) filter = 1, or comment vanilla ddp and cbf filter-related commands
% 5) plots and tables: comment what is not needed 
% 6) number of DDP iterations (will impact computational load)

clear all; close all; clc;
tStart = tic;

% Add Paths
parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory)); 

%% Initialize system parameters

% Numerical parameters
dt = 0.02;  % Sampling
N = 150;    % Horizon (# of iteration that discretizes each pass)
T = 0:dt:dt*N-1*dt;

% Set final state
xf = [1;1];

% Set range of initial values to accept and generate a set of x0
size_of_set = 50;

%% Safety Constraints
[safety_fun,obs_loc,circ] = obstacles_2d();

x0 = [];
while length(x0) < size_of_set
    lim = 0.6*xf(1);
    ic = lim*rand(2,1) - 0; % Generates random points in the first 25%
    % of the nominal path
    collision = (obs_loc(1) - ic(1))^2 + (obs_loc(2) - ic(2))^2;  % Obstacle
    if collision > 2.25*(obs_loc(3)^2)   % Double-checks if the IC is in or 
        % too close to the obstacle
% % %         if ic(1) < 0.9*ic(2) || ic(1) > 1.1*ic(2) % Makes sure x!=y for nominal 
% % % %       trajectory, if obstacle is centered at x-y
        x0 = [x0, ic];
    end
end

% Quick breakpoint to check the location of ICs
figure()
rectangle('Position',circ(1,:),'Curvature',[1 1],'LineWidth', 2, 'LineStyle','--'); hold on; grid on;
for kk = 1:length(x0)
    plot(x0(1,kk),x0(2,kk),'*','LineWidth',2)
end
axis equal; axis square;

%% Barrier Conditions
% % % % % % % % alpha_h = 10;
% % % alpha_h = 20; up to 25 % Best tune for CBF-Aware, backwards pass only
% % % alpha_h = 35; % Tuned value for CBF outside DDP

h = safety_fun.h;
hx = safety_fun.hx;
% % % % % % % % % alpha_cbf = @(x) alpha_h*h{1}(x);
% % % % % % % % % filter = 1; % Code is not ready to run without the filter. If needed, 
% % % % % % % % % % comment the lines of the propagation of states, checks of collision and
% target, and table of comparison.
%% Call System Dynamics

f_dyn = single_integrator_dynamics(dt,1);
n = length(f_dyn.x); % State dimensions
m = length(f_dyn.u); % Input dimensions
f_dyn.x_des = xf; f_dyn.u_des = [0;0];

%% Quadratic costs (running cost and terminal cost)

Q = 0.05*eye(n);    % States weights
R = 1e-2*eye(m);    % Controls weights
S = 10*eye(n);      % Final cost

% % % % Best tune for CBF-Aware, backwards pass only @ alpha = 20 to 25
% % % Q = 0.05*eye(n);    % States weights --- Works decent up to 0.15*eye
% % % R = 1e-2*eye(m);    % Controls weights
% % % S = 10*eye(n);      % Final cost

% Cost functions
run_cost = @(x,u,deriv_bool) run_quad_cost(x, u, Q, R, xf, deriv_bool);
term_cost = @(x,deriv_bool) terminal_quad_cost(x, xf, S, deriv_bool);

%% Nominal input and state

xbar = []; ubar = [];
% Takes initial condition and nominal control and propagates the states 
% according to the system's dynamics
for i = 1:length(x0)
    ubar = 0.0*ones(m, N-1,i);  % nominal control
    xbar(:,1,i) = x0(:,i);      % initial state
    for k = 1:N-1
        xbar(:,k + 1,i) = f_dyn.F(xbar(:,k,i), ubar(:,k,i));
    end
end

%% Optimization parameters

iter = 100;             % Maximum number of iterations
toler = 1e-3;           % For cost change

% For the regularization part by Yuichiro Aoyama
lambda = 0.5;            % Initial value for lambda for regularization
dlambda= 1;              % Initial value for dlambda for regularization
lambdaFactor = 1.6;      % Lambda scaling factor
lambdaMax = 1e10;        % Lambda maximum value
lambdaMin = 1e-6;        % Below this value lambda = 0

opt_par = struct('iter',iter,'toler',toler,'lambda', lambda, 'dlambda',...
    dlambda,'lambdaFactor', lambdaFactor, 'lambdaMax', lambdaMax,...
    'lambdaMin', lambdaMin);

ddp_2nd_order = 0; % For Full DDP

%% Discrete DDP with Safety Penalties (CBF-Aware DDP)

for i = 1:length(x0)

    f_dyn.x0 = x0(:,i);
    sys_par = struct('dt', dt, 'N', N,'x0', x0(:,i), 'xf', xf, 'n', n, 'm', m);
    
    % Propagate states over time with CBF-Aware DDP
    [X(:,:,i), U(:,:,i), ~, ~, ~, ~, ~, ~, ~, ii(i), ~, J(i), safety(i), w_back(:,i), w_fwd(:,i)] = safe_disc_ddp_alg(ddp_2nd_order,...
        f_dyn, run_cost, term_cost, sys_par, ubar(:,:,i), xbar(:,:,i), opt_par, h, hx, obs_loc);

    % Propagate states with vanilla DDP + CBF filter outside for comparison
    [X_ddp(:,:,i),U_ddp(:,:,i), ~, ~,~, ~, ~, K_u, ii_ddp(i), ~, ~] = disc_ddp_alg(ddp_2nd_order,...
        f_dyn, run_cost, term_cost, sys_par, ubar(:,:,i), xbar(:,:,i), opt_par);
% % % % % % %     [Xcbf(:,:,i), Ucbf(:,:,i), Jcbf(i),safety2(i)] = CBF_Filter(X_ddp(:,:,i), U_ddp(:,:,i), K_u, f_dyn, run_cost,...
% % % % % % %         term_cost, sys_par, hx, alpha_cbf, obs_loc);

    % Check if the particles reached the target
    if X(1,end,i) >= 0.99*xf(1) && X(1,end,i) <= 1.025*xf(1) && X(2,end,i) >= 0.99*xf(2) && X(2,end,i) <= 1.025*xf(2)
        target(i) = 1;
    else
        target(i) = 0;
    end
% % % % % % %     if Xcbf(1,end,i) >= 0.99*xf(1) && Xcbf(2,end,i) >= 0.99*xf(2)
% % % % % % %         target_cbf(i) = 1;
% % % % % % %     else
% % % % % % %         target_cbf(i) = 0;
% % % % % % %     end
end

%% Plottings
% Paths
figure()
rectangle('Position',circ(1,:),'Curvature',[1 1],'LineWidth', 2, 'LineStyle','--');
hold on; grid on;
invisible = line(NaN, NaN, 'LineWidth', 2, 'LineStyle', '--', 'Color',[0 0 0]); % Workaround to have obstacle's info displayed
plot(xf(1),xf(2),'*','LineWidth',2)
title('X-Y path (top view)','Interpreter','latex');
ylabel('$y$','FontName','Times New Roman', 'Interpreter','latex');
xlabel('$x$','FontName','Times New Roman', 'Interpreter','latex');
leg = ['Center at (' num2str(obs_loc(1)) ', ' num2str(obs_loc(2)) ') and radius = ' num2str(obs_loc(3))];
legend(leg,'Location','northwest','FontName','Times New Roman', 'Interpreter','latex');
for i = 1:length(x0)
% % %     plot(x0(1,i),x0(2,i),'*','LineWidth',2)
    plot(X(1,:,i),X(2,:,i),'HandleVisibility','off')
% % % % % % %     plot(Xcbf(1,:,i),Xcbf(2,:,i),'--')
end
box on; axis square; axis equal; 
% xlim = [-1 1]; ylim = [-1 1];

% Controls
% U1 (CBF-Aware)
figure()
hold on; grid on;
title('Control vs time','Interpreter','latex');
ylabel('$u_{1}$','FontName','Times New Roman', 'Interpreter','latex');
xlabel('$t$ in seconds','FontName','Times New Roman', 'Interpreter','latex');
for i = 1:size_of_set
    plot(T(1:end-1), U(1,:,i));
end

% U2 (CBF-Aware)
figure()
hold on; grid on;
title('Control vs time','Interpreter','latex');
ylabel('$u_{1}$','FontName','Times New Roman', 'Interpreter','latex');
xlabel('$t$ in seconds','FontName','Times New Roman', 'Interpreter','latex');
for i = 1:size_of_set
    plot(T(1:end-1), U(2,:,i));
end

% w/alpha (CBF-Aware)
% Backwards
figure()
hold on; grid on;
title('$\alpha$ from the backward propagation vs time','Interpreter','latex');
ylabel('$\alpha$ backwards','FontName','Times New Roman', 'Interpreter','latex');
xlabel('$t$ in seconds','FontName','Times New Roman', 'Interpreter','latex');
for i = 1:size_of_set
    plot(T(1:end), w_back(:,i));
end

% Forwards
figure()
hold on; grid on;
title('$\alpha$ from the forward propagation vs time','Interpreter','latex');
ylabel('$\alpha$ forwards','FontName','Times New Roman', 'Interpreter','latex');
xlabel('$t$ in seconds','FontName','Times New Roman', 'Interpreter','latex');
for i = 1:size_of_set
    plot(T(1:end-1), w_fwd(:,i));
end

%% Check just a few controls
% % % % % figure()
% % % % % hold on; grid on;
% % % % % title('Control vs time','FontName','Times New Roman','Interpreter','latex');
% % % % % ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
% % % % % xlabel('Time (s)','FontName','Times New Roman','Interpreter','latex');
% % % % % plot(T(1:end-1), U(1,:,231),'b');
% % % % % plot(T(1:end-1), U(2,:,231),'r');
% % % % % plot(T(1:end-1), Ucbf(1,:,231),'--b');
% % % % % plot(T(1:end-1), Ucbf(2,:,231),'--r');
% % % % % plot(T(1:end-1), U_ddp(1,:,231),':b');
% % % % % plot(T(1:end-1), U_ddp(2,:,231),':r');
% % % % % legend('$u_{1}$ CBF-Aware DDP', '$u_{2}$ CBF-Aware DDP', '$u_{1}$ CBF Filtering',...
% % % % %     '$u_{2}$ CBF Filtering', '$u_{1}$ vanilla','$u_{2}$ vanilla',...
% % % % %     'FontName', 'Times New Roman', 'Interpreter','latex');

%% Performance comparison

IC = []; Cost = []; Iterations = []; Collision = [];
Target = []; Cost2 = []; Iterations2 = []; Collision2 = []; Target2 = [];
NoCollisionCost = []; NoCollisionCost2 = []; NoCollision_IC = [];
NoCollisionTarget = []; NoCollisionTarget2 = [];

for i = 1:length(x0)
    IC = [IC; num2str(i, '%03d')];
    Cost = [Cost; J(i)]; 
% % % % %     Cost2 = [Cost2; Jcbf(i)];saveas(figure(1), [parentDirectory '/cbf-aware_ddp/results/5p_optdecay.png']);
    Iterations = [Iterations; ii(i)];
% % % % %     Iterations2 = [Iterations2; ii_ddp(i)];
    Collision = [Collision; safety(i)];
% % % % %     Collision2 = [Collision2; safety2(i)];
    Target = [Target; target(i)];
% % % % %     Target2 = [Target2; target_cbf(i)];
% % % % %     % Filter out ICs that did not collide
% % % % %     if safety(i) == 0 && safety2(i) == 0
% % % % %         NoCollision_IC = [NoCollision_IC; num2str(i, '%03d')];
% % % % %         NoCollisionCost = [NoCollisionCost; J(i)];
% % % % %         NoCollisionCost2 = [NoCollisionCost2; Jcbf(i)];
% % % % %         NoCollisionTarget = [NoCollisionTarget; target(i)];
% % % % %         NoCollisionTarget2 = [NoCollisionTarget2; target_cbf(i)];
% % % % %     end
end

compare = table(IC, Cost, Iterations, Collision, Target);
% % % % % compareCBFoutside = table(IC, Cost2, Iterations2, Collision2, Target2);
% % % % % T1 = outerjoin(compare,compareCBFoutside,'MergeKeys', true);
% % % % % T2 = mergevars(T1,{'Cost','Iterations','Collision','Target'},...
% % % % %                'NewVariableName','CBF-Aware DDP','MergeAsTable',true);
% % % % % T3 = mergevars(T2,{'Cost2','Iterations2','Collision2','Target2'},...
% % % % %                'NewVariableName','CBF as a Filter (Outside DDP)','MergeAsTable',true);
% % % % % disp(T3);
disp(compare);

mean_cost = mean(Cost); % mean_cost2 = mean(Cost2);
total_collisions = sum(Collision); % total_collisions2 = sum(Collision2);
total_target = sum(Target); % total_target2 = sum(Target2);

% TO-DO: 1) Write headers and write a row with # of collisions, # of target
% reached and mean cost; 2) Get the clean results inside the table; 
% 3) Add timer???

% Enter a name that grabs params and saves figures, tables and data
% consistently --- TO-DO: automate that
writetable(compare,[parentDirectory '/cbf-aware_ddp/results/' num2str(size_of_set) 'p_optdecay-4.xlsx']);
save([parentDirectory '/cbf-aware_ddp/results/' num2str(size_of_set) 'p_optdecay-4.mat']);
saveas(figure(1), [parentDirectory '/cbf-aware_ddp/results/' num2str(size_of_set) 'p_optdecay_path-4.png']);
saveas(figure(2), [parentDirectory '/cbf-aware_ddp/results/' num2str(size_of_set) 'p_optdecay_u1-4.png']);
saveas(figure(3), [parentDirectory '/cbf-aware_ddp/results/' num2str(size_of_set) 'p_optdecay_u2-4.png']);
saveas(figure(4), [parentDirectory '/cbf-aware_ddp/results/' num2str(size_of_set) 'p_optdecay_alphaback-4.png']);
saveas(figure(5), [parentDirectory '/cbf-aware_ddp/results/' num2str(size_of_set) 'p_optdecay_alphafwd-4.png']);

tEnd = toc(tStart);
fprintf('Total elapsed time (min): %.2f\n',tEnd/60);
fprintf('Obstacle centered at (%.2f, %.2f) with a radius of %.2f\n',obs_loc(1), obs_loc(2), obs_loc(3));
fprintf('\nMean Cost of CBF-Aware: %.4f\n',mean_cost);
% % % % % fprintf('Mean Cost of CBF Filter: %.4f\n',mean_cost2);
% % % % % fprintf('Cost improvement of %.2f%%\n',((mean_cost2-mean_cost)/mean_cost2)*100); % DIVIDES BY COST OR COST2?
fprintf('Number of safety violations using CBF-Aware: %d (%.2f%%)\n',total_collisions,(total_collisions/size_of_set)*100);
% % % % % fprintf('Number of safety violations using CBF Filter: %d (%.2f%%)\n',total_collisions2,(total_collisions2/size_of_set)*100);
fprintf('Number of targets reached using CBF-Aware: %d (%.2f%%)\n',total_target,(total_target/size_of_set)*100);
% % % % % fprintf('Number of targets reached using CBF Filter: %d (%.2f%%)\n',total_target2,(total_target2/size_of_set)*100);
% % % % % 
% % % % % 
% % % % % no_collision_compare = table(NoCollision_IC, NoCollisionCost, NoCollisionCost2, NoCollisionTarget, NoCollisionTarget2);
% % % % % mean_NoCollisionCost = mean(NoCollisionCost); mean_NoCollisionCost2 = mean(NoCollisionCost2);
% % % % % debug_target = sum(NoCollisionTarget); debug_target2 = sum(NoCollisionTarget2);
% % % % % 
% % % % % writetable(no_collision_compare,[parentDirectory '/cbf-aware_ddp/results/250p_a10_q005_xf4-4_AssymetricalObstacle-NOCOLLISION.xlsx']);
% % % % % 
% % % % % fprintf('\nMean Collision-free Cost of CBF-Aware: %.4f\n',mean_NoCollisionCost);
% % % % % fprintf('Mean Collision-free Cost of CBF Filter: %.4f\n',mean_NoCollisionCost2);
% % % % % fprintf('Cost improvement of %.2f%%\n',((mean_NoCollisionCost2-mean_NoCollisionCost)/mean_NoCollisionCost2)*100);
% % % % % fprintf('\nNumber of Collision-free targets reached using CBF-Aware: %d\n',debug_target);
% % % % % fprintf('Number of Collision-free targets reached using CBF Filter: %d\n',debug_target2);