function rho = solve_rho(u, W, x, theta, eps, d_theta, N_theta) %% rho与x有关，与theta无关 
 %% 插值准备
     d_theta_f = 1e-4;
     [X, Theta] = meshgrid(x, theta);
     theta_f = 0:d_theta_f:1;
     [X_f, Theta_f] = meshgrid(x, theta_f);
     %% spline插值
     u_aux = u;
     W_aux = [W, W(:, 1)];
     u_f = interp1(theta, u_aux, theta_f, 'spline');
     W_f = interp2(X, Theta, W_aux', X_f, Theta_f, 'spline');
     W_f = W_f';  
     rho = (theta_f(2)-theta_f(1)) * sum(W_f(:, 1:end-1) .* exp(u_f(1:end-1)/eps), 2);

end

