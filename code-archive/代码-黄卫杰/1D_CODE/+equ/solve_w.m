function W = solve_w(u, W, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta)

pL_u = equ.WENO5_left_RXR(u, d_theta);
pR_u = equ.WENO5_right_RXR(u, d_theta);

du_m_dw = 0 .* W;
for kk = 1: size(W, 1)
     pL_W = equ.WENO5_left_RXR(W(kk, :), d_theta);
     pR_W = equ.WENO5_right_RXR(W(kk, :), d_theta);
     du_m_dw(kk, :) = equ.RF_flux_neg_uv(pR_u, pL_u, pR_W, pL_W);
end

%%%系数矩阵
A_diag2 = diag(K - rho, 0);
H_mat = reshape(H, [1, 1, N_theta]);
A_H_cell = squeeze(num2cell(H_mat.*eye(N_x+1), [1, 2]));
A_H = blkdiag(A_H_cell{:});
A = Ap - dt * kron(diag(ones(N_theta, 1), 0), A_diag2) +...
    - dt * A_H;

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

% 右端项
% b = eps * W - dt * 2 * eps * du_m_dw;
b = eps * W - dt * b_part1 - dt * b_part2; 
b = reshape(b, [], 1);
%向前Euler    
W_new = A \ b;
W = reshape(W_new, size(W));% 更新
end

