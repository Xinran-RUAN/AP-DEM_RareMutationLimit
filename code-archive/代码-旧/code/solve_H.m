<<<<<<< HEAD
function [H] = solve_H(N_theta, N_x, D, Tri_C, K, rho)
%% 系数矩阵C表示，D(theta),
%% 一个theta对应一个D(theta)，对应1个特征值，所以需要循环
H = zeros(1, N_theta);   
for kk = 1: N_theta % 计算特征值
    C = - D(kk) .* Tri_C - diag(K(2:N_x) - rho(2: N_x));
    [~, H(kk)] = eigs(C, 1, 'smallestreal');
end

end

=======
function [H] = solve_H(N_theta, N_x, D, Tri_C, K, rho)
%% 系数矩阵C表示，D(theta),
%% 一个theta对应一个D(theta)，对应1个特征值，所以需要循环
H = zeros(1, N_theta);   
for kk = 1: N_theta % 计算特征值
    C = - D(kk) .* Tri_C - diag(K(2:N_x) - rho(2: N_x));
    [~, H(kk)] = eigs(C, 1, 'smallestreal');
end

end

>>>>>>> d725b20c4e0dc455060136e6e0d3a79bba8f525a
