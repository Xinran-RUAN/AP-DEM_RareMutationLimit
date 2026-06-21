function [u_inte, W_inte] = uW_spline_theta(u, W, x, theta, theta_f)
[X, Theta] = meshgrid(x, theta);
u = u;
W = [W, W(:, 1)];
[X_f, Theta_f] = meshgrid(x, theta_f);
 
u_inte = interp1(theta, u, theta_f, 'spline');
W_inte = interp2(X, Theta, W', X_f, Theta_f, 'spline');
W_inte = W_inte';
end

