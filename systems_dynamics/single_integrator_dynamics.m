%% Dynamics of Planar Double Integrator
function f_dyn = single_integrator_dynamics(dt,deriv_bool)
%     syms x1 x2 u1 u2
%     x=[x1;x2];
%     u=[u1;u2];
%     n=2; m=2;
    n = 2; m = 2; % state and input dimentions
    x = sym('x',[n 1]);
    u = sym('u',[m 1]);
% x1dot= u1;
% x2dot= u2;
% xdot=f(x,u)
% f=[u(1); u(2)];
    
    % Discrete dynamics
    f_dyn.F=@(x,u) [x(1);x(2)] + dt*[u(1); u(2)];
    
    f_dyn.f = [0; 0]; % Should I multiply by dt???
    f_dyn.g = eye(2); % Should I multiply by dt???

%     [fx, fu, fxx, fxu, fuu]=deal([]);

    if deriv_bool
        % dynamics gradients
%         fx=[0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
        f_dyn.fx=@(x,u)  eye(n)+dt*[0 0;0 0];
%         fu=[0 0;0 0;1 0;0 1];
        f_dyn.fu=@(x,u)  dt*[1 0;0 1];
% % %         fx=@(x,u) dt*[0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
% % %         fu=@(x,u) dt*[0 0;0 0;1 0;0 1];

        fx_temp = f_dyn.fx(x,u);
        for ii=1:n
                fxx_temp = jacobian(fx_temp(ii,:), x);
                f_dyn.fxx{ii}=matlabFunction(fxx_temp,'Vars',{x,u});
                fxu_temp = jacobian(fx_temp(ii,:), u).';
                f_dyn.fxu{ii}=matlabFunction(fxu_temp,'Vars',{x,u});
        end
        fu_temp = f_dyn.fu(x,u);
        for ii=1:m
                fuu_temp = jacobian(fu_temp(ii,:), u);
                f_dyn.fuu{ii}=matlabFunction(fuu_temp,'Vars',{x,u});
        end
    end
    
    f_dyn.x=x; f_dyn.u=u; f_dyn.n=n; f_dyn.m=m;

end


