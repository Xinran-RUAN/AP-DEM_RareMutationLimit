function info = mass_balance_stats_1d(rho, Kx, op)
%MASS_BALANCE_STATS_1D 计算一维质量平衡所需的空间积分量。
%
% 对模型 eps*dM/dt = int rho*(K-rho) dx，定义：
%   M   = int rho dx
%   R   = int rho*(K-rho) dx
%   A   = int K*rho dx
%   B   = int rho^2 dx
%
% rho 是空间总密度 rho(x)，不是 trait marginal P(theta)。

rho = real(rho(:));
Kx = real(Kx(:));
if numel(Kx) ~= numel(rho)
    error('mass_balance_stats_1d:SizeMismatch', 'rho 和 Kx 长度不一致。');
end
wx = utils.trapz_weights(op.nx, op.dx);
wx = wx(:);
if numel(wx) ~= numel(rho)
    error('mass_balance_stats_1d:WeightMismatch', '积分权重长度与 rho 长度不一致。');
end

info = struct();
info.mass = real(wx.' * rho);
info.reactionIntegral = real(wx.' * (rho .* (Kx - rho)));
info.A = real(wx.' * (Kx .* rho));
info.B = real(wx.' * (rho.^2));
info.rhoMin = min(rho);
info.rhoMax = max(rho);
end
