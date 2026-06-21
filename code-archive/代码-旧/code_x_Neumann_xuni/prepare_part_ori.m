function [A_con, A_pre] = prepare_part_ori(eps, dx, dt, dtheta, D, Nx, Nth)

%% 对x的二阶导部分，关于x是Neumann边界条件
beta = dt * eps^2 / dtheta^2;
alpha = dt .* D ./ dx^2;
alpha = reshape(alpha, 1, 1, []);
A_pre = diag(ones(Nx+1, 1), 0) .* (eps+ 2 * alpha + 2 * beta) + ...
    - diag(ones(Nx, 1), 1) .* alpha + ...
    - diag(ones(Nx, 1), -1) .* alpha;
A_pre(1, 2, :) = - 2 * alpha;
A_pre(Nx+1, Nx, :) = - 2 * alpha;

%% 对theta的二阶导部分，关于theta是周期边界条件
A_sub = -beta .* diag(ones(Nx+1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_con = blkdiag(A_cell{:});
A_con = A_con + kron(diag(ones(Nth-1, 1), 1), A_sub) ...
              + kron(diag(ones(Nth-1, 1), -1), A_sub);
A_con(end-Nx:end, 1:Nx+1) = A_sub;
A_con(1:Nx+1, end-Nx:end) = A_sub;

end

