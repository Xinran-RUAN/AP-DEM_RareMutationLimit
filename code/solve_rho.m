function rho = solve_rho(u, W, x, theta, eps, d_theta, N_theta) %% rho与x有关，与theta无关 
 %% 插值准备
d_theta_f = 1e-3;
[X, Theta] = meshgrid(x, [theta, 1]);
theta_f = 0:d_theta_f:1; 
[X_f, Theta_f] = meshgrid(x, theta_f);
%% spline插值  
u_aux = [u, u(:, 1)];  
w_aux = [W, W(:, 1)];  
uinte = interp1([theta, 1], u_aux, theta_f, 'spline');
winte = interp2(X, Theta, w_aux', X_f, Theta_f, 'spline');
winte = winte';
rho = (theta_f(2)-theta_f(1)) * sum(winte(:, 1:end-1) .* exp(uinte(1:end-1)/eps), 2);

% rho = d_theta * sum(W(:, 1:N_theta) .* exp(u(1:N_theta)/eps), 2);

end

