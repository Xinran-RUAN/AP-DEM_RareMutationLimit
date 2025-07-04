function W = solve_w(u, W, W_temp, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta, theta, ol)
u_plus = [u(2:end), u(1)];
u_minus = [u(end), u(1:end-1)];

if ol == 0
    pL = (u - u_minus) / d_theta; % 左差分 u_j - u_{j-1}
    pR = (u_plus - u) / d_theta; % 右差分 u_{j+1} - u_j
    id1 = (pL .* pR >= 0) .* (abs(pL) <= abs(pR));
    id2 = (pL .* pR >= 0) .* (abs(pL) > abs(pR));
    id3 = pL .* pR < 0;
    du = pL .* id1 + pR .* id2 + 0 .* id3;
    du(1) = pR(1);
     % du = weno5_diff(u, d_theta);
else
    if ol == 1
        %% spline插值求导数
        du = 0 .* u;
        pp = spline(theta(2:end), u(2:end));
        pp_der = fnder(pp, 1);
        du(2:end) = ppval(pp_der, theta(2:end));
        pR = (u_plus - u) / d_theta; % 右差分 u_{j+1} - u_j
        du(1) = pR(1);
    else
        pp = spline([theta, 1], [u, u(1)]);
        pp_der = fnder(pp, 1);
        du = ppval(pp_der, theta);
    end
end
    
if ol == 0
    %% 迎风
    W_theta_plus = W_temp(:, [2:end, 1]); % 关于theta是周期边界条件
    W_theta_minus = W_temp(:, [end, 1:end-1]);
    WL = (W_temp - W_theta_minus) / d_theta;
    WR = (W_theta_plus - W_temp) / d_theta; % 右差分，dudtheta<0，
    
    dw = (du <= 0) .* WL + (du > 0) .* WR;
    dw(:, 1) = WR(:, 1);
%     dw = 0 .* W;
%     for ii = 1: size(W, 1)
%         dw(ii, :) = weno5_diff(W(ii, :), d_theta);
%     end
    
else
    if ol == 1
        dw = 0 .* W;
        %% spline 插值求导数
        for ii = 1: size(W, 1)
            pp = spline(theta(2:end), W(ii, 2:end));
            pp_der = fnder(pp, 1);
            dw(ii, 2:end) = ppval(pp_der, theta(2:end));
        end
        W_theta_plus = W_temp(:, [2:end, 1]);
        WR = (W_theta_plus - W_temp) / d_theta; % 右差分，dudtheta<0，
        dw(:, 1) = WR(:, 1);
    else
        dw = 0 .* W;
        for ii = 1: size(W, 1)
            pp = spline([theta, 1], [W(ii, :), W(ii, 1)]);
            pp_der = fnder(pp, 1);
            dw(ii, :) = ppval(pp_der, theta);
        end
    end
end

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

