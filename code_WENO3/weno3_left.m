function phix_minus = weno3_left(phi, dx)
% WENO3 差分型负向导数重构
N = length(phi);
phix_minus = zeros(size(phi));
epsilon = 1e-1;

for i = 1:N
    i_m2 = mod(i-3, N) + 1;
    i_m1 = mod(i-2, N) + 1;
    i_p0 = i;
    i_p1 = mod(i, N) + 1;

    % 差分
    a = (phi(i_p0) - phi(i_m1)) / dx; % φ_i - φ_{i-1}
    b = (phi(i_m1) - phi(i_m2)) / dx; % φ_{i-1} - φ_{i-2}
    c = (phi(i_p1) - phi(i_p0)) / dx; % φ_{i+1} - φ_i

    % 子模板导数
    phi0 = (3/2)*a - (1/2)*b;  % stencil {i-2,i-1,i}
    phi1 = (1/2)*a + (1/2)*c;  % stencil {i-1,i,i+1}

    % 光滑性指标
    IS0 = (phi(i_p0) - 2*phi(i_m1) + phi(i_m2))^2;
    IS1 = (phi(i_p1) - 2*phi(i_p0) + phi(i_m1))^2;

    % 权重
    alpha0 = (1/3) / (epsilon + IS0)^2;
    alpha1 = (2/3) / (epsilon + IS1)^2;
    alpha_sum = alpha0 + alpha1;
    w0 = alpha0 / alpha_sum;
    w1 = alpha1 / alpha_sum;

    % 合成导数
    phix_minus(i) = w0 * phi0 + w1 * phi1;
end
end
