function W = solve_w(u, W, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta)

pL_u = WENO5_left_RXR(u, d_theta);
pR_u = WENO5_right_RXR(u, d_theta);

du_m_dw = 0 .* W;
for kk = 1: size(W, 1)
     pL_W = WENO5_left_RXR(W(kk, :), d_theta);
     pR_W = WENO5_right_RXR(W(kk, :), d_theta);
     du_m_dw(kk, :) = RF_flux_neg_uv(pR_u, pL_u, pR_W, pL_W);
end

%%%系数矩阵
A_diag2 = diag(K - rho, 0);
H_mat = reshape(H, [1, 1, N_theta]);
A_H_cell = squeeze(num2cell(H_mat.*eye(N_x+1), [1, 2]));
A_H = blkdiag(A_H_cell{:});
A = Ap - dt * kron(diag(ones(N_theta, 1), 0), A_diag2) +...
    - dt * A_H;

% 右端项
b = eps * W - dt * 2 * eps * du_m_dw;
b = reshape(b, [], 1);
%向前Euler    
W_new = A \ b;
W = reshape(W_new, size(W));% 更新
end

