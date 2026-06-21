function [H] = solve_H_fem(N_theta, mesh_divide_elements, mesh_divide_nodes, D, K, rho, xi, eta, weight, K1, K2_x, K2_y)
%% 系数矩阵C表示，D(theta),
%% 一个theta对应一个D(theta)，对应1个特征值，所以需要循环
H = zeros(1, N_theta); 

for kk = 1: N_theta % 计算特征值
    D_k = D(kk);
    [A, B] = stiffness_matrix_H(mesh_divide_elements, mesh_divide_nodes, D_k, K, rho, xi, eta, weight, K1, K2_x, K2_y);
    [~, H(kk)] = eigs(A, B, 1, 'smallestreal');
end  

end

