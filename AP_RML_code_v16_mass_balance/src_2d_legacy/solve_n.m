function ne = solve_n(epsilon, dt, dx, dtheta, Nx, Nth, ne, D, K)
%% 
beta = epsilon^2 * dt / dtheta^2;
alpha = dt .* D ./ dx^2;
alpha = reshape(alpha, 1, 1, []);
A_pre = diag(ones(Nx-1, 1), 0) .* (1+ 2 * alpha + 2 * beta) + ...
    - diag(ones(Nx-2, 1), 1) .* alpha + ...
    - diag(ones(Nx-2, 1), -1) .* alpha;
A_pre(1, 1, :) = 1 + alpha + 2 * beta;
A_pre(end, end, :) = 1 + alpha + 2 * beta;
A_sub = -beta .* diag(ones(Nx-1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_con = blkdiag(A_cell{:});
A_con = A_con + kron(diag(ones(Nth-1, 1), 1), A_sub) + kron(diag(ones(Nth-1, 1), -1), A_sub);
A_con(end-Nx+2:end, 1:Nx-1) = A_sub;
A_con(1:Nx-1, end-Nx+2:end) = A_sub;

rho = dtheta * sum(ne(:, 1:Nth), 2);
    %% 线性插值
    % ne_aux = [ne, ne(:, 1)];
    % neinte = interp2(xx, thth, ne_aux', xxinte, ththinte, 'spline');
    % neinte = neinte';
    % rho = 0.001 * sum(neinte(:, 1:end-1), 2);
 
    %% solve the equation to obtain n(x_j,\theta_i,t_m+1)
    %%%系数矩阵   
    A_diag2 = diag(K(2:Nx) - rho(2:Nx), 0); 
    A = A_con - dt * kron(diag(ones(Nth, 1), 0), A_diag2);

    % 右端项  
    b = ne; 
    b = reshape(b(2: Nx, :), [], 1);
    %向前Euler       
    ne_new = A \ b;               
    ne = reshape(ne_new, size(ne(2:Nx, :)));% 更新
    ne = [ne(1, :); ne; ne(end, :)];
    
    
end

