function f_dyn_bar=Safety_Embedding_dynamics(f_dyn,bas_dyn)
    n=f_dyn.n;  x=f_dyn.x;   u=f_dyn.u;  m =f_dyn.m; nbas=bas_dyn.nbas;
    u_des=f_dyn.u_des;
    % augment bas dynamics
    f_dyn_bar.x = sym('x',[n+nbas 1]); %n_bas is number of barrier states;
    f_dyn_bar.F=@(x,u)[f_dyn.F(x,u);bas_dyn.f_z(x,u)];
    f_dyn_bar.f=@(x)[f_dyn.f(x);bas_dyn.f_z(x)];
    f_dyn_bar.g=@(x)[f_dyn.g(x);bas_dyn.g_z(x)];
    
    xbar0=[f_dyn.x0;bas_dyn.z0]; xbar_des=[f_dyn.x_des;bas_dyn.zf];
    nbar=length(xbar0);

    f_dyn_bar.fx=@(x,u)[f_dyn.fx(x,u) zeros(n,nbas);bas_dyn.fx_z(x,u)];
    f_dyn_bar.fu=@(x,u)[f_dyn.fu(x,u);bas_dyn.fu_z(x,u)];

    f_dyn_bar.A=f_dyn_bar.fx(xbar_des,u_des);
    f_dyn_bar.B=f_dyn_bar.fu(xbar_des,u_des);

    f_dyn_bar.x0 = xbar0; f_dyn_bar.x_des = xbar_des;
    f_dyn_bar.n = nbar; f_dyn_bar.m = m;
end