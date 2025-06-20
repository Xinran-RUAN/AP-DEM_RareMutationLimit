function [H] = solve_H(N_theta, N_x, D, Tri_C, K, rho)
H = zeros(1, N_theta);
opt.tol = 1e-14;
opt.maxit = 1e4;
for kk = 1: N_theta % ��������ֵ
    C = - D(kk) .* Tri_C - diag(K(2:N_x) - rho(2: N_x));
    [~, H(kk)] = eigs(C, 1, 'smallestreal', opt);
%     if ~(all(VV>=0) || all(VV<=0))
%         VV;
%     end
end
end

