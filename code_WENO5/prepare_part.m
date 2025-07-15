function [B, Ap, Tri_C] = prepare_part(eps, dt, dx, d_theta, N_x, N_theta, D)

beta = eps * dt / d_theta^2;  
B = (1 + 2 * beta) * diag(ones(N_theta, 1), 0) +...
            - beta * diag(ones(N_theta-1, 1), 1) +...
            - beta * diag(ones(N_theta-1, 1), -1);
B(1, N_theta) = - beta;
B(N_theta, 1) = - beta;  

alpha_x = dt .* D ./ dx^2;
alpha_theta = eps^2 * dt / d_theta^2;
alpha_x = reshape(alpha_x, 1, 1, []); % alpha为三阶张量，beta为标量
A_pre = diag(ones(N_x-1, 1), 0) .* (eps+2*alpha_x+2*alpha_theta) + ...
      - diag(ones(N_x-2, 1), 1) .* alpha_x + ...
      - diag(ones(N_x-2, 1), -1) .* alpha_x;  
A_pre(1, 1, :) = eps + alpha_x + 2 * alpha_theta;
A_pre(end, end, :) = eps + alpha_x + 2 * alpha_theta;
A_sub = - alpha_theta .* diag(ones(N_x - 1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_diag = blkdiag(A_cell{:});
Ap = A_diag + kron(diag(ones(N_theta - 1, 1), 1), A_sub) +...
            + kron(diag(ones(N_theta - 1, 1), -1), A_sub);
Ap(end-N_x+2:end, 1:N_x-1) = A_sub;
Ap(1:N_x-1, end-N_x+2:end) = A_sub;

gamma = 1 / dx^2;
Tri_C = - 2 * gamma * diag(ones(N_x - 1, 1), 0) +...
            + gamma * diag(ones(N_x - 2, 1), 1) +...
            + gamma * diag(ones(N_x - 2, 1), -1);
Tri_C(1, 1) = - gamma;
Tri_C(end, end) = - gamma;
end

