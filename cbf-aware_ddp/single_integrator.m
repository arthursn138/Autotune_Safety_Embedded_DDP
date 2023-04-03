%% Control-Limited Differential Dynamic Programming
% Arthur Nascimento, Hassan Almubarak
% nascimento@gatech.edu, halmubarak@gatech.edu
% ACDS Lab @ Georgia Tech
% Last Update June 28 2022

% This is an example for 2D single integrator.
% This example uses optimal decay to allow for a better choice of alpha.

% TO-DO: UPDATE INSTRUCTIONS BELOW

% Before running, make sure the obstacle is where you want it to be, and
% your alpha_h and Q are well tuned. Lower alpha_h will avoid the obstacle
% before (more conservative); Lower Q will pull the safe trajectory more
% to the nominal, i.e. low Q avoids later, but can diverge more easily.
% Good values for the obstacle at a favorable location are alpha_h in
% between 10 and 25; Q in between 0.05*eye and 0.1*eye. (For CBF-Aware).
% For CBF filter, you have to go a little higher on alpha_h (~35)

clear all; clc;
% % % clear all; close all; clc;

% % Add Paths
% general_path = '../Autotune_Safety_Embedded_DDP';
% addpath(genpath(general_path))  
parentDirectory = fileparts(cd);
addpath(genpath(parentDirectory)); 

%% Initialize system parameters
t_total = tic;
% Numerical parameters
dt = 0.02;  % Sampling
N = 150;    % Horizon (# of iteration that discretizes each pass)
T = 0:dt:dt*N-1*dt;
% T = dt*N;

% Set initial and final state
x0 = [0;0];
xf = [3;2.5];

% System Dynamics
f_dyn = single_integrator_dynamics(dt,1);

n = length(f_dyn.x); % State dimensions
m = length(f_dyn.u); % Input dimensions

f_dyn.x0 = x0; f_dyn.x_des = xf; f_dyn.u_des = [0;0]; 

sys_par=struct('dt', dt, 'N', N,'x0', x0, 'xf', xf, 'n', n, 'm', m);

%% Safety Constraints
[safety_fun,obs_loc,circ] = obstacles_2d();

% Button to toggle/untoggle the simulation of CBF outside DDP ("CBF-filered
% DDP") for comparison purposes. 1 or True = on; 0 or False = off.
filter = 1;

%% Barrier Conditions
alpha_h = 7;
% % % alpha_h = 20; % Best tune for CBF-Aware, backwards pass only
% % % alpha_h = 35; % Tuned value for CBF outside DDP
% % % Higher aplha makes it more conservative, i.e. dodge before and higher
% % % cost

h = safety_fun.h;
hx = safety_fun.hx;
alpha_cbf = @(x) alpha_h*h{1}(x);

%% Quadratic costs (running cost and terminal cost)

Q = 0.05*eye(n);    % States weights
R = 1e-2*eye(m);    % Controls weights
S = 10*eye(n);      % Final cost

% % % Q = 0.05*eye(n);    % States weights ---  Best tune for CBF-Aware, backwards pass only
% % % R = 1e-2*eye(m);    % Controls weights ---  Best tune for CBF-Aware, backwards pass only
% % % S = 10*eye(n);      % Final cost ---  Best tune for CBF-Aware, backwards pass only
% % % 
% % % Those values work great for x0 = (0,0); xf = (1,1); obstacle r = 0.1,
% % % x_center = 0.4, y_center = 0.5!!!!!!!!!!!!!!!

% Cost functions
run_cost = @(x,u,deriv_bool) run_quad_cost(x, u, Q, R, xf, deriv_bool);
term_cost = @(x,deriv_bool) terminal_quad_cost(x, xf, S, deriv_bool);

%% Nominal input and state

ubar = 0.0*ones(m, N-1);    % nominal control
xbar = []; 
xbar(:,1) = x0;             % initial state

% Takes initial condition and nominal control and propagates the states 
% according to the system's dynamics
for k = 1:N-1
    xbar(:,k + 1) = f_dyn.F(xbar(:,k), ubar(:,k));
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

%% Discrete Vanilla DDP
tic
[X_ddp,U_ddp, ~, ~,~, ~, ~, K_u, ii_ddp, ~, J_ddp] = disc_ddp_alg(ddp_2nd_order,...
    f_dyn, run_cost, term_cost, sys_par, ubar, xbar, opt_par);

tf_ddp = toc;

%% CBF Filtering
if filter
    [Xcbf, Ucbf, Jcbf, safety_cbf] = CBF_Filter(X_ddp, U_ddp, K_u, f_dyn, run_cost,...
        term_cost, sys_par, hx, alpha_cbf, obs_loc);
    tf_ddp = toc;
else
    Xcbf = []; Ucbf = []; Jcbf = [];
    tf_ddp = toc;
end

%% Discrete DDP with Safety Penalties (CBF-Aware DDP)
tic

[X, U, ~, ~, alphaLinsearch, ~, ~, ~, ~, ii, ~, J, safety, w_back, w_fwd] = safe_disc_ddp_alg(ddp_2nd_order,...
    f_dyn, run_cost, term_cost, sys_par, ubar, xbar, opt_par, h, hx, obs_loc);

tf = toc;
fprintf('Number of DDP iterations: %d \n', ii);

%% Plottings

if filter == 0
    Xcbf = []; Ucbf = [];
end
X_ddp = []; U_ddp = [];

si_cbf_plot(X, Xcbf, X_ddp, U, Ucbf, U_ddp, T, x0, xf, circ, filter, w_back, w_fwd, obs_loc);
% % % saveas(figure(1), [parentDirectory '/cbf-aware_ddp/results/ReplicatingManan.png']);

% % % % % saveas(figure(1), [parentDirectory '/cbf-aware_ddp/results/MinNorm_CompareFilter_DiffCoeff-Path.png']);
% % % % % saveas(figure(2), [parentDirectory '/cbf-aware_ddp/results/MinNorm_CompareFilter_DiffCoeff-U.png']);
% % % % % saveas(figure(3), [parentDirectory '/cbf-aware_ddp/results/MinNorm_CompareFilter_DiffCoeff-AlphaBack.png']);
% % % % % saveas(figure(4), [parentDirectory '/cbf-aware_ddp/results/MinNorm_CompareFilter_DiffCoeff-AlphaFwd.png']);

tf_total = toc(t_total);
fprintf('Total runtime: %.2f minutes\n', tf_total/60);
fprintf('Total cost: %.4f\n', J);
%% Performance comparison
% % % % % % 
% % % % % % % % % % Find min and max disturbance for all the controls
% % % % % % % % % [minU1, maxU1] = find_min_max(U(1,:), 0.08);
% % % % % % % % % [minU2, maxU2] = find_min_max(U(2,:), 0.047);
% % % % % % % % % [minU1cbf, maxU1cbf] = find_min_max(Ucbf(1,:), 0.08);
% % % % % % % % % [minU2cbf, maxU2cbf] = find_min_max(Ucbf(2,:), 0.05);
% % % % % % % % % 
% % % % % % % % % maxU1diff = (maxU1cbf-maxU1);
% % % % % % % % % maxU2diff = (maxU2cbf-maxU2);
% % % % % % % % % minU1diff = (minU1-minU1cbf);
% % % % % % % % % minU2diff = (minU2-minU2cbf);
% % % % % % 
% Generate table with results
if filter
    Method = ["CBF-Aware DDP"; "CBF filtered DDP"];
    FinalCost = [J; Jcbf];
    FinalState1 = [X(1,end); Xcbf(1,end)];
    FinalState2 = [X(2,end); Xcbf(2,end)];
else
    Method = ["CBF-Aware DDP"; "Vanilla DDP"];
    FinalCost = [J; J_ddp];
    FinalState1 = [X(1,end); X_ddp(1,end)];
    FinalState2 = [X(2,end); X_ddp(2,end)];
end

Iterations = [ii; ii_ddp];
Time = [tf; tf_ddp];
compare = table(Method, FinalCost, Iterations, Time, FinalState1, FinalState2);

% Display results
disp(compare);
fprintf('Cost improvement of %.2f%%\n',((Jcbf-J)/Jcbf)*100); % Is the division correct?
% % % % % % % % % % % % % % fprintf('Max control effort difference of %.2f for U1, and %.2f for U2\n', maxU1diff, maxU2diff);
% % % % % % % % % % % % % % fprintf('Max control effort improvement of %.2f%% for U1, %.2f%% for U2\n', 100*maxU1diff/maxU1cbf, 100*maxU2diff/maxU2cbf);
% % % % % % % % % % % % % % fprintf('Min control effort difference of %.2f for U1, %.2f for U2\n', minU1diff, minU2diff);
% % % % % % % % % % % % % % fprintf('Min control effort improvement of %.2f%% for U1, %.2f%% for U2\n', 100*maxU1diff/maxU1cbf, 100*maxU2diff/maxU2cbf);
% % % % % % % % % % % 
% % % % % % % fprintf('\nQP in Backwards + FWD propagations\n');
% % % % % % % % % % % % % % % fprintf('Objective functions: Back = u*Quu*u + Qu*u; FWD = minimum norm\n');
% % % % % % % fprintf('Objective functions: Back = minimum norm; FWD = minimum norm\n');
% % % % % % % fprintf('\nParameters:\n');
% % % % % % % fprintf('alpha_h = %f\n', alpha_h);
% % % % % % % fprintf('Q = %.2f*I\n', Q(1));
% % % % % % % fprintf('R = %.2f\n', R(1));
% % % % % % % fprintf('Obstacle at (%.2f, %.2f) with radius = %.2f\n', obs_loc(1), obs_loc(2), obs_loc(3));
% % % % % % % fprintf('Initial condition = (%.2f, %.2f), Target = (%.2f, %.2f)\n', x0(1), x0(2), xf(1),xf(2));
% % % % % % % % % % % 

% % % % % % % % % save([parentDirectory '/cbf-aware_ddp/results/ReplicatingManan.mat']);

%% Auxiliary functions
% % % % % % 
% % % % % % function [min_v, max_v] = find_min_max(v, threshold)
% % % % % % % Calculate and compare control effort
% % % % % % % Isolate for when spikes due to obstacle occurs and get indexes
% % % % % % % Vector "v" should be a collumn or row vector, matrices are not accepted
% % % % % % % Threshold is a percentage and should be written in the decimal form from
% % % % % % % 0 to 1., i.e. 5% should be written as 0.05
% % % % % % 
% % % % % % growth_rate = zeros(1,length(v)); index = []; j = 0;
% % % % % % for j = 2:length(v)
% % % % % %     growth_rate(j) = abs(1-v(j-1)/v(j));
% % % % % %     if growth_rate(j) > threshold
% % % % % %         index= [index, j];
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % % Find min and max disturbance. (robust to the optimal control being a
% % % % % % % constant ascent or descent)
% % % % % % minv =[]; maxv = []; regime = [];
% % % % % % for k = index(1)-1:index(end)+1 % Give some room to look for what is off-optimal
% % % % % %     rate = v(k)-v(k-1);
% % % % % %     if rate < 0 % descent
% % % % % %         regime(k) = 0;
% % % % % %     elseif rate == 0 % flat
% % % % % %         regime(k) = 1;
% % % % % %     else  % ascent
% % % % % %         regime(k) = 2;
% % % % % %     end
% % % % % %     % Lists all changes of descent/ascent/flat
% % % % % %     if regime(k-1) < regime(k)
% % % % % %         minv = [minv; v(k-1)];
% % % % % %     elseif regime(k-1) > regime(k)
% % % % % %         maxv = [maxv; v(k-1)];
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % min_v = min(minv);
% % % % % % max_v = max(maxv);
% % % % % % 
% % % % % % end
