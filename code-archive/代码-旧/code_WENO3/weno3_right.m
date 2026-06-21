function phix_plus = weno3_right(phi, dx)
% WENO3 差分型正向导数重构
% phi: 1D 数组，网格点的函数值
% dx: 网格间距
% phix_plus: 正向导数的 WENO3 近似

N = length(phi);
phix_plus = zeros(size(phi));
eps0 = 1e-1;

for i = 1:N
    % 周期索引
    i_m1 = mod(i-2, N) + 1; % i-1
    i_p0 = i;               % i
    i_p1 = mod(i, N) + 1;   % i+1
    i_p2 = mod(i+1, N) + 1; % i+2

    % 差分
    a = (phi(i_p1) - phi(i_p0)) / dx; % Δ?φ_i
    b = (phi(i_p0) - phi(i_m1)) / dx; % Δ?φ_{i-1}
    c = (phi(i_p2) - phi(i_p1)) / dx; % Δ?φ_{i+1}

    % 两个子模板导数近似
    phi0 = (1/2)*a + (1/2)*b;
    phi1 = (3/2)*a - (1/2)*c;

    % 光滑性指标
    IS0 = (a - b)^2;
    IS1 = (a - c)^2;

    % 权重（线性权重：phi0 对应 2/3，phi1 对应 1/3）
    alpha0 = (2/3) / (eps0 + IS0)^2;
    alpha1 = (1/3) / (eps0 + IS1)^2;
    alpha_sum = alpha0 + alpha1;
    w0 = alpha0 / alpha_sum;
    w1 = alpha1 / alpha_sum;

    % 合成导数
    phix_plus(i) = w0 * phi0 + w1 * phi1;
end
end
