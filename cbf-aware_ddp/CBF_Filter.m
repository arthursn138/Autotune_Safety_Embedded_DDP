function [xCBF, uCBF, Jcbf, safetycbf] = CBF_Filter(x_star, u_star, K_u, f_dyn,...
    run_cost, term_cost, ddp_par, hx, alpha_cbf, obs_loc)
% CBF filter using quadratic programming
% Set QPsolver as 1 for fmincon or 0 for quadprog

% Initializations:
N = ddp_par.N; n = ddp_par.n; m = ddp_par.m; x0 = ddp_par.x0;
xCBF = zeros(n,N); uCBF = zeros(m,N-1); L = 0; safetycbf = zeros(1,N);
xCBF(:,1) = x0; K_u = flip(K_u,3);
fmincon_options = optimoptions(@fmincon,'Display','off');

% Main Loop

for i = 1:N-1
    deltax = xCBF(:,i) - x_star(:,i);
    uCBF(:,i) = u_star(:,i) + K_u(:,:,i)*deltax;
    ku = uCBF(:,i);

    hx_cbf = hx{1}(x_star(1:2,i));

    fminimize = @(uCBF) (uCBF - ku)'*(uCBF - ku);
    uCBF(:,i) = fmincon(fminimize, ku, [], [], [], [], [], [],...
        @(uCBF) constraint_cbf(uCBF, hx_cbf, f_dyn.f, f_dyn.g, alpha_cbf,...
        x_star(:,i)), fmincon_options);

    xCBF(:,i+1) = f_dyn.F(xCBF(:,i), uCBF(:,i));
    L = L + run_cost(xCBF(:,i),uCBF(:,i),false);
    safetycbf(i + 1) = (obs_loc(1) - xCBF(1,i+1))^2 + (obs_loc(2) - xCBF(2,i+1))^2;
end

Jcbf = L + term_cost(xCBF(:,end),true);
if min(safetycbf(2:end)) <= obs_loc(3)^2
    safetycbf = 1;
else
    safetycbf = 0;
end
end


function [c, ceq] = constraint_cbf(ku,hx_cbf,f,g,alpha_cbf,x)
c = -(hx_cbf*(f + g*ku) + alpha_cbf(x));
ceq = [];
end