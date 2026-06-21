function [A1, A2, A3, A4,...
          b1, b2] = element_stiffness_matrix_W(K1, k2_x, k2_y,...
                               Con, x_local, y_local, xi, eta, weight, D, H,...
                               W_local, dumdw, K_local, rho_local,...
                                      eps, dt, dtheta)

Jacobi = (x_local(2) - x_local(1)) * (y_local(3) - y_local(1)) - (x_local(3) - x_local(1)) * (y_local(2) - y_local(1));
xi_x = (y_local(3) - y_local(1)) / Jacobi;
xi_y = (x_local(1) - x_local(3)) / Jacobi;

eta_x = (y_local(1) - y_local(2)) / Jacobi;
eta_y = (x_local(2) - x_local(1)) / Jacobi;

A1_point = eps * K1 .* Jacobi;
A2_point = dt * D * (k2_x + k2_y) * Jacobi;
A3_point = - dt * eps^2 * K1 * Jacobi / dtheta^2;

K_point = K_local(1) * varphi_1_2D(xi, eta, 1) +...
        + K_local(2) * varphi_1_2D(xi, eta, 2) + ...
        + K_local(3) * varphi_1_2D(xi, eta, 3);
rho_point = rho_local(1) * varphi_1_2D(xi, eta, 1) +...
          + rho_local(2) * varphi_1_2D(xi, eta, 2) + ...
          + rho_local(3) * varphi_1_2D(xi, eta, 3);

K_point_e = permute(repmat(K_point, 3, 1, 3), [3, 1, 2]);
rho_point_e = permute(repmat(rho_point, 3, 1, 3), [3, 1, 2]);
 
A4_point = - dt * (H * K1 + K1.*(K_point_e - rho_point_e)) .* Jacobi;

W_point = W_local(1) * varphi_1_2D(xi, eta, 1) +...
        + W_local(2) * varphi_1_2D(xi, eta, 2) + ...
        + W_local(3) * varphi_1_2D(xi, eta, 3);
W_point_e = permute(repmat(W_point, 1, 1, 3), [3, 1, 2]);
b1_point = eps * Con .* W_point_e * Jacobi;
dumdw_point = dumdw(1) * varphi_1_2D(xi, eta, 1) +...
            + dumdw(2) * varphi_1_2D(xi, eta, 2) + ...
            + dumdw(3) * varphi_1_2D(xi, eta, 3);
dumdw_point_e = permute(repmat(dumdw_point, 1, 1, 3), [3, 1, 2]);
b2_point = - dt * dumdw_point_e .* Con * Jacobi;

weight_point = permute(repmat(weight, 3, 1, 3), [3, 1, 2]);
A1 = sum(A1_point .* weight_point, 3);
A2 = sum(A2_point .* weight_point, 3);
A3 = sum(A3_point .* weight_point, 3);
A4 = sum(A4_point .* weight_point, 3);

weight_point2 = permute(repmat(weight, 1, 1, 3), [3, 1, 2]);
b1 = sum(b1_point .* weight_point2, 3);
b2 = sum(b2_point .* weight_point2, 3);
end

