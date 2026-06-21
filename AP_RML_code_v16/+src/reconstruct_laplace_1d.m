function [rho, n, info] = reconstruct_laplace_1d(W, u, eps, op, allowFallback)
%RECONSTRUCT_LAPLACE_1D  周期三点 Laplace 重构 rho。
%
%  该函数是 WKB 框架的一部分：仍然使用 W 和 u 来近似
%      rho(x) = int W(x,theta) exp(u(theta)/eps) dtheta。
%  它并不求解原始密度方程。
%
%  小 epsilon 时必须避免直接计算 Wstar*exp(Ustar/eps) 造成溢出。
%  因此这里在 log 变量中计算：
%      log rho_j = log Wstar_j + Ustar/eps
%                  + 0.5*log(2*pi*eps/|u_{theta theta}(theta*)|).
%
%  若二次峰值不是凹峰，或插值得到的 Wstar 非正，则回退到 direct-log。
%  回退也只是 WKB 重构，不是 direct density solver。

if nargin < 5
    allowFallback = true;
end

peak = src.phase_peak_info_1d(u, op);
info = struct();
info.mode = 'laplace';
info.usedFallback = false;
info.thetaStar = peak.thetaStar;
info.dduStar = peak.dduStar;
info.uStar = peak.uStar;
info.uNodeMax = peak.uNodeMax;
info.uStarMinusNodeMax = peak.uStar - peak.uNodeMax;
info.logrhoMax = NaN;
info.numBadWstar = 0;
info.numCappedHigh = 0;

if ~peak.ok || ~isfinite(peak.dduStar) || peak.dduStar >= -1e-14
    if allowFallback
        [rho, n, infoFallback] = src.reconstruct_rho_1d(W, u, eps, op, 'direct-log');
        info = local_merge(info, infoFallback);
        info.usedFallback = true;
        info.mode = 'direct-log-fallback-bad-curvature';
        return;
    else
        error('Laplace reconstruction failed: nonnegative or nearly zero second derivative at the maximum.');
    end
end

km = peak.km; k0 = peak.k0; kp = peak.kp;
L = peak.L;
Wstar = real(W(:,km))*L(1) + real(W(:,k0))*L(2) + real(W(:,kp))*L(3);
info.numBadWstar = nnz(~isfinite(Wstar) | Wstar <= 0);

if any(~isfinite(Wstar)) || any(Wstar <= 0)
    if allowFallback
        [rho, n, infoFallback] = src.reconstruct_rho_1d(W, u, eps, op, 'direct-log');
        info = local_merge(info, infoFallback);
        info.usedFallback = true;
        info.mode = 'direct-log-fallback-bad-Wstar';
        return;
    else
        error('Laplace reconstruction failed: interpolated W at the phase maximum is non-positive.');
    end
end

logPref = 0.5*log(2*pi*eps/abs(peak.dduStar));
logrho = log(max(Wstar, realmin)) + peak.uStar/eps + logPref;
logCap = 690;  % exp(690) 仍在 double 安全范围内。
info.logrhoMax = max(logrho(:));
info.logrhoMin = min(logrho(:));
info.numCappedHigh = nnz(logrho > logCap);

rho = exp(min(logrho, logCap));
rho = real(rho(:));

% n 只用于图像和 gauge-invariant 诊断。这里仍用 direct-log 的稳定实现。
[~, n] = src.reconstruct_rho_1d(W, u, eps, op, 'direct-log');
end

function out = local_merge(a,b)
out = a;
if ~isstruct(b), return; end
f = fieldnames(b);
for i = 1:numel(f)
    out.(f{i}) = b.(f{i});
end
end
