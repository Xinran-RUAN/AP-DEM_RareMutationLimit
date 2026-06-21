function rho = solve_rho(u, W, theta, eps, d_theta, N_theta, ol) %% rho与x有关，与theta无关 
 %% 插值准备
 if ol == 0
     dtheta_f = 1e-4;
     theta_f = 0:dtheta_f:1;
     theta_ext = [theta, 1];
     u_ext = [u, u(1)];
     u_f = interp1(theta_ext, u_ext, theta_f, 'spline');
     
     %% W的spline插值
     W_ext = [W, W(:, 1)];
     Wf_tmp = interp1(theta_ext, W_ext', theta_f, 'spline');
     W_f = Wf_tmp';  
     rho = (theta_f(2)-theta_f(1)) * sum(W_f(:, 1:end-1) .* exp(u_f(1:end-1)/eps), 2);
 else
     rho = d_theta * sum(W(:, 1:N_theta) .* exp(u(1:N_theta)/eps), 2);
 end
end

