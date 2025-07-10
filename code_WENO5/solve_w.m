function W = solve_w(u, W, W_temp, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta, theta, ol)

 pL_u = WENO5_left_RXR(u, d_theta);
 pR_u = WENO5_right_RXR(u, d_theta);

du_m_dw = 0 .* W;

for kk = 1: size(W, 1)
     pL_W = WENO5_left_RXR(W(kk, :), d_theta);
     pR_W = WENO5_right_RXR(W(kk, :), d_theta);
    du_m_dw(kk, :) = RF_flux_neg_uv(pR_u, pL_u, pR_W, pL_W);
end


%%%系数矩阵
A_diag2 = diag(K(2:N_x) - rho(2:N_x), 0);
H_mat = reshape(H, [1, 1, N_theta]);
A_H_cell = squeeze(num2cell(H_mat.*eye(N_x-1), [1, 2]));
A_H = blkdiag(A_H_cell{:});
A = Ap - dt * kron(diag(ones(N_theta, 1), 0), A_diag2) +...
    - dt * A_H;

% 右端项
b = eps * W - dt * 2 * eps * du_m_dw;
b = reshape(b(2: N_x, :), [], 1);
%向前Euler  
w_new = A \ b; 
W = reshape(w_new, size(W(2:N_x, :)));% 更新
W = [W(1, :); W; W(end, :)];
    
end

