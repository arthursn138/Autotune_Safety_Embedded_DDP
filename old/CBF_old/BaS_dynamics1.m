%% Dynamics of Continous Time Barrier States (BaS)
function bas_dyn=BaS_dynamics1(f_dyn,safety_fun,gamma,dt)
    x=f_dyn.x; u=f_dyn.u; n=f_dyn.n; 
    x0=f_dyn.x0; xf=f_dyn.x_des; uf=f_dyn.u_des; 
    F=f_dyn.F(x,u); f=f_dyn.f(x,u);  g=f_dyn.g(x,u); 
    h=safety_fun.h; hx=safety_fun.hx; 
    hxx=safety_fun.hxx; nbas=safety_fun.nbas;

    % define BaS and its initial and final states
    beta=0;
    for ii=1:length(h)
        beta=beta+1/h{ii}(x);
    end
    beta = matlabFunction(beta,'Vars',{x});
    beta_0=beta(xf);
    z0=beta(x0) - beta_0;
    zf=beta(xf) - beta_0;
    
    x = sym('x',[n+nbas 1]); %n_bas is number of barrier states;
    z = sym('z',[nbas 1]); %n_bas is number of barrier states;
    
    z=x(n+nbas)+beta_0;
    
    phi_0 = @(zeta) - zeta^2;
    phi_1 = @(zeta,eta) eta*zeta^2 - zeta;
    
    for ii=1:length(h)
    zeta = z - (beta(f_dyn.x)-1/h{ii}(f_dyn.x));
    eta = h{ii}(f_dyn.x);
    F_z = phi_0(zeta)*hx{ii}(f_dyn.x)*F - gamma*phi_1(zeta,eta);
    f_z = phi_0(zeta)*hx{ii}(f_dyn.x)*f - gamma*phi_1(zeta,eta);
    g_z = phi_0(zeta)*hx{ii}(f_dyn.x)*g;
    end

    fx_z=jacobian(F_z,x);
    fu_z=jacobian(F_z,u); 

    % symbolic to function handle
    F_z= matlabFunction(F_z,'Vars',{x,u});
    fx_z= matlabFunction(fx_z,'Vars',{x,u});
    fu_z= matlabFunction(fu_z,'Vars',{x,u});
    f_z= matlabFunction(f_z,'Vars',{x,u});
    g_z= matlabFunction(g_z,'Vars',{x,u});
    
    bas_dyn.F_z=F_z;
    bas_dyn.f_z=f_z;
    bas_dyn.g_z=g_z;
    bas_dyn.fx_z=fx_z;
    bas_dyn.fu_z=fu_z;
    bas_dyn.z0=z0;
    bas_dyn.zf=zf;
    bas_dyn.nbas=nbas;
end
