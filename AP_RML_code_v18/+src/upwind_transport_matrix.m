function T = upwind_transport_matrix(u, op)
%UPWIND_TRANSPORT_MATRIX Matrix for D_theta F, F=-W u_theta.
% Vector ordering is [all x nodes at theta_1; all x nodes at theta_2; ...].
%
% 向量化实现：先构造 theta 方向 KxK 稀疏矩阵，再 kron(Ix)。
% 旧实现逐块 append rows/cols/vals，profile 中会产生明显开销。

K = op.Ntheta;
nx = op.nx;
dth = op.dtheta;
u = real(u(:).');

rows = zeros(3*K,1);
cols = zeros(3*K,1);
vals = zeros(3*K,1);
p = 0;

for k = 1:K
    kp = mod(k, K) + 1;
    km = mod(k-2, K) + 1;

    a_kph = -(u(kp) - u(k)) / dth;
    a_kmh = -(u(k) - u(km)) / dth;

    c0 = max(a_kph,0)/dth - min(a_kmh,0)/dth;
    cp = min(a_kph,0)/dth;
    cm = -max(a_kmh,0)/dth;

    p = p + 1; rows(p) = k; cols(p) = k;  vals(p) = c0;
    p = p + 1; rows(p) = k; cols(p) = kp; vals(p) = cp;
    p = p + 1; rows(p) = k; cols(p) = km; vals(p) = cm;
end

Sth = sparse(rows(1:p), cols(1:p), vals(1:p), K, K);
T = kron(Sth, op.Ix);
end
