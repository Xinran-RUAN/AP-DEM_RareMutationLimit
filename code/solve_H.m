<<<<<<< HEAD
function [H] = solve_H(N_theta, N_x, D, Tri_C, K, rho)
%% ϵ������C��ʾ��D(theta),
%% һ��theta��Ӧһ��D(theta)����Ӧ1������ֵ��������Ҫѭ��
H = zeros(1, N_theta);   
for kk = 1: N_theta % ��������ֵ
    C = - D(kk) .* Tri_C - diag(K(2:N_x) - rho(2: N_x));
    [~, H(kk)] = eigs(C, 1, 'smallestreal');
end

end

=======
function [H] = solve_H(N_theta, N_x, D, Tri_C, K, rho)
%% ϵ������C��ʾ��D(theta),
%% һ��theta��Ӧһ��D(theta)����Ӧ1������ֵ��������Ҫѭ��
H = zeros(1, N_theta);   
for kk = 1: N_theta % ��������ֵ
    C = - D(kk) .* Tri_C - diag(K(2:N_x) - rho(2: N_x));
    [~, H(kk)] = eigs(C, 1, 'smallestreal');
end

end

>>>>>>> d725b20c4e0dc455060136e6e0d3a79bba8f525a
