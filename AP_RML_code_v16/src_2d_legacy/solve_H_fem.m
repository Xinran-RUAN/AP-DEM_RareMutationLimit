function [H] = solve_H_fem(N_theta, mesh_divide_elements, mesh_divide_nodes, D, K, rho, xi, eta, weight, K1, K2_x, K2_y)
%SOLVE_H_FEM Principal generalized eigenvalue for each theta node.
%   Uses the smallest algebraic eigenvalue of A v = H B v.

H = zeros(1, N_theta);
opts = struct();
opts.tol = 1e-11;
opts.maxit = 500;

for kk = 1:N_theta
    D_k = D(kk);
    [A, B] = stiffness_matrix_H(mesh_divide_elements, mesh_divide_nodes, D_k, K, rho, xi, eta, weight, K1, K2_x, K2_y);
    A = (A + A')/2;
    B = (B + B')/2;
    try
        [~, lam] = eigs(A, B, 1, 'sa', opts);
    catch
        try
            [~, lam] = eigs(A, B, 1, 'smallestreal', opts);
        catch
            lamAll = eig(full(A), full(B));
            lam = min(real(lamAll));
        end
    end
    H(kk) = real(lam);
end
end
