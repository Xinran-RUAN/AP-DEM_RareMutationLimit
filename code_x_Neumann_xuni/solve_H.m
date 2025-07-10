function [H] = solve_H(N_theta, N_x, D, Tri_C, K, rho)
%% 系数矩阵C表示，D(theta),
%% 一个theta对应一个D(theta)，对应1个特征值，所以需要循环
H = zeros(1, N_theta);   
for kk = 1: N_theta % 计算特征值
    C = - D(kk) .* Tri_C - diag(K - rho);
    [~, H(kk)] = eigs(C, 1, 'smallestreal');
end  

end

