function rho = spline_integration(eps, theta, u, W)
     d_theta_f = 1e-4;
     theta_f = 0:d_theta_f:1;
     %% spline²åÖµ
     u_f = interp1(theta, u', theta_f, 'spline');
     W_f = interp1(theta, W', theta_f, 'spline');
     rho = d_theta_f * sum(W_f(:, 1:end-1) .* exp(u_f(1:end-1)/eps), 2);

end

