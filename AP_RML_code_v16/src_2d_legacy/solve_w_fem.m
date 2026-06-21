function W = solve_w_fem(u, W, K, rho, D, H, eps, d_theta, dt, N_theta,...
                        mesh_divide_elements, mesh_divide_nodes,...
                        xi, eta, weight, K1, K2_x, K2_y, Con)
                    
pL_u = WENO5_left_RXR(u, d_theta);
pR_u = WENO5_right_RXR(u, d_theta);

du_m_dw = 0 .* W;
for kk = 1: size(W, 1)
     pL_W = WENO5_left_RXR(W(kk, :), d_theta);
     pR_W = WENO5_right_RXR(W(kk, :), d_theta);
     du_m_dw(kk, :) = RF_flux_neg_uv(pR_u, pL_u, pR_W, pL_W);
end
du_m_dw = du_m_dw * 2 * eps;

%% partial_theta W partial u的拆分，以下是右端项中的2*eps*partial_theta(-W \partial_theta u)
up = [u(2:end), u(1)];
um = [u(end), u(1:end-1)];
a_kph = -(up - u)/d_theta; % a_{k+1/2}
a_kmh = -(u - um)/d_theta; % a_{k-1/2}
Wp = [W(:, 2:end), W(:, 1)];
Wm = [W(:, end), W(:, 1:end-1)];
F_hat_kph = max(a_kph, 0) .* W + min(a_kph, 0) .* Wp;
F_hat_kmh = max(a_kmh, 0) .* Wm + min(a_kmh, 0) .* W;

b_part1 = 2 * eps * (F_hat_kph - F_hat_kmh) / d_theta;

%% partial_theta W partial u的拆分，以下是右端项2eps*W*ddu，还未移到右端
up = [u(2:end), u(1)];
um = [u(end), u(1:end-1)];
ddu = (up - 2*u + um) / d_theta^2;
b_part2 = 2 * eps * W .* ddu; 

du_m_dw = b_part1 + b_part2;

[A, b] = stiffness_matrix_W(K1, K2_x, K2_y, Con, mesh_divide_elements,...
                            mesh_divide_nodes, xi, eta, weight,...
                            D, W, du_m_dw, K, rho, H, eps, dt, d_theta, N_theta);

%向前Euler      
W_new = A \ b;
W = reshape(W_new, size(W));% 更新
end

