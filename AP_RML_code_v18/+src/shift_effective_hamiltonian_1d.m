function [Huse, info] = shift_effective_hamiltonian_1d(Hraw, par)
%SHIFT_EFFECTIVE_HAMILTONIAN_1D  对有效 Hamiltonian 做标量 gauge 平移。
%
%  WKB 分解中，H 只以如下成对方式出现：
%       u_t + \hat H - eps*u_{theta theta} = -H,
%       eps*W_t + ... = W*(... + H + ...).
%  因此把 H 替换为 H-C(t)，并同时在相位和振幅方程中使用同一个
%  H-C(t)，不会改变重构密度 n=W*exp(u/eps) 的连续方程。这个 C(t)
%  只是 WKB 常数 gauge 的选择。
%
%  小 epsilon 下，如果直接使用带有较大公共偏移的 H，则 h-only 积分因子
%  exp(dt*H/eps) 会把这个公共偏移也放大，导致 W 出现无意义的整体指数
%  缩放。默认选择 C=min(H)，使 min(Huse)=0。这样 fittest trait 对应
%  零惩罚，其它 trait 的相位相对下降，同时显著减少 eps 刚性。
%
%  注意：这里只做“标量”平移，不做 theta 依赖的 rephase，因此不会改变
%  相位方程的导数结构，也不会人为改变 thetaWKB。

if nargin < 2 || isempty(par)
    par = struct();
end
Hraw = real(Hraw(:).');
mode = 'min';
if isfield(par, 'hamiltonianGauge') && ~isempty(par.hamiltonianGauge)
    mode = lower(char(par.hamiltonianGauge));
end

switch mode
    case {'min','minimum','zero-min'}
        C = min(Hraw);
    case {'mean','average'}
        C = mean(Hraw);
    case {'max','maximum'}
        C = max(Hraw);
    case {'none','off','raw'}
        C = 0;
    otherwise
        error('Unknown hamiltonianGauge "%s".', mode);
end

Huse = Hraw - C;
info = struct();
info.hamiltonianGaugeMode = mode;
info.hamiltonianGaugeShift = C;
info.HrawMin = min(Hraw);
info.HrawMax = max(Hraw);
info.HrawRange = max(Hraw)-min(Hraw);
info.HuseMin = min(Huse);
info.HuseMax = max(Huse);
info.HuseRange = max(Huse)-min(Huse);
end
