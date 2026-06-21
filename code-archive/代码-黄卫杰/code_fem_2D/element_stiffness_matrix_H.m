function [A_local, B_local] = element_stiffness_matrix_H(k1, k2_x, k2_y, x_local, y_local, K_local, rho_local, xi, eta, weight, D)
Jacobi = (x_local(2) - x_local(1)) * (y_local(3) - y_local(1)) - (x_local(3) - x_local(1)) * (y_local(2) - y_local(1));

xi_x = (y_local(3) - y_local(1)) / Jacobi;
xi_y = (x_local(1) - x_local(3)) / Jacobi;

eta_x = (y_local(1) - y_local(2)) / Jacobi;
eta_y = (x_local(2) - x_local(1)) / Jacobi;

K_point = K_local(1) * varphi_1_2D(xi, eta, 1) +...
        + K_local(2) * varphi_1_2D(xi, eta, 2) + ...
        + K_local(3) * varphi_1_2D(xi, eta, 3);
rho_point = rho_local(1) * varphi_1_2D(xi, eta, 1) +...
          + rho_local(2) * varphi_1_2D(xi, eta, 2) + ...
          + rho_local(3) * varphi_1_2D(xi, eta, 3);

K_point_e = permute(repmat(K_point, 3, 1, 3), [3, 1, 2]);
rho_point_e = permute(repmat(rho_point, 3, 1, 3), [3, 1, 2]);

A_local_point = (D * (k2_x + k2_y) - k1.*(K_point_e - rho_point_e)) .* Jacobi;
B_local_point = k1 .* Jacobi;

weight_point = permute(repmat(weight, 3, 1, 3), [3, 1, 2]);
A_local = sum(A_local_point .* weight_point, 3);
B_local = sum(B_local_point .* weight_point, 3);
end

