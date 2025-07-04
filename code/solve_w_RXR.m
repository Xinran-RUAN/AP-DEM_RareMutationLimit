function W = solve_w_RXR(u, W, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta)
du_plus = WENO5_right(u, d_theta);
pR = reshape(du_plus, size(u));
du_minus = WENO5_left(u, d_theta);
pL = reshape(du_minus, size(u));

id1 = (pL .* pR >= 0) .* (abs(pL) <= abs(pR));
id2 = (pL .* pR >= 0) .* (abs(pL) > abs(pR));
id3 = pL .* pR < 0;
du = pL .* id1 + pR .* id2 + 0 .* id3;


    
%% 迎风
W_theta_plus = W(:, [2:end, 1]); % 关于theta是周期边界条件
W_theta_minus = W(:, [end, 1:end-1]);
WL = (W - W_theta_minus) / d_theta;
WR = (W_theta_plus - W) / d_theta; % 右差分，dudtheta<0，

dw = (du <= 0) .* WL + (du > 0) .* WR;
dw(:, 1) = WR(:, 1);
    

%%%系数矩阵
A_diag2 = diag(K(2:N_x) - rho(2:N_x), 0);
H_mat = reshape(H, [1, 1, N_theta]);
A_H_cell = squeeze(num2cell(H_mat.*eye(N_x-1), [1, 2]));
A_H = blkdiag(A_H_cell{:});
A = Ap - dt * kron(diag(ones(N_theta, 1), 0), A_diag2) +...
    - dt * A_H;

% 右端项
b = eps * W + dt * 2 * eps * dw .* du;
b = reshape(b(2: N_x, :), [], 1);
%向前Euler
w_new = A \ b;
W = reshape(w_new, size(W(2:N_x, :)));% 更新
W = [W(1, :); W; W(end, :)];
    
end

