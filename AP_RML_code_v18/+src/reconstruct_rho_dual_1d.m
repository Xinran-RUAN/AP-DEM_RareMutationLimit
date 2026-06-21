function [rho, nU, info, WU] = reconstruct_rho_dual_1d(Wc, uU, eps, opW, opU, mode, interpMethod)
%RECONSTRUCT_RHO_DUAL_1D 双 theta 网格 WKB 重构。
% Wc: nx x Kw，在粗 theta_W 网格；uU: 1 x Ku，在细 theta_u 网格。
% 先把 W 插值到 u 细网格，再做 log-stabilized quadrature。
if nargin < 6 || isempty(mode), mode = 'direct-log'; end
if nargin < 7 || isempty(interpMethod), interpMethod = 'pchip'; end

info = struct('mode', mode, 'dualThetaGrid', true, 'WtoUInterp', interpMethod, ...
              'numNonpositiveW', 0, 'lognMax', NaN, 'logrhoMax', NaN, ...
              'numCappedNHigh', 0, 'numCappedRhoHigh', 0, 'usedFallback', false);

% 目前双网格推进默认用 direct-log。若用户写 laplace/hybrid，仍用细 u 网格
% 上的 log quadrature；后处理中再做 Laplace/phase curvature 重构。
modeLower = lower(char(mode));
if any(strcmp(modeLower, {'laplace','laplace-hybrid','hybrid-laplace'}))
    info.usedFallback = true;
    info.requestedMode = mode;
    modeLower = 'direct-log';
end

switch modeLower
    case {'direct','direct-log'}
        WU = src.interp_periodic_theta_1d(opW.theta, Wc, opU.theta, interpMethod);
        Wreal = real(WU);
        bad = (~isfinite(Wreal)) | (Wreal <= 0);
        info.numNonpositiveW = nnz(bad);
        Wpos = max(Wreal, realmin);
        logW = log(Wpos);
        logn = bsxfun(@plus, logW, real(uU(:).')/eps);
        info.lognMax = max(logn(:));
        info.lognMin = min(logn(:));
        info.numNonfiniteLogN = nnz(~isfinite(logn(:)));
        logNCap = 690;
        info.numCappedNHigh = nnz(logn > logNCap);
        nU = exp(min(logn, logNCap));
        nU(~isfinite(nU)) = 0;
        weights = opU.dtheta * ones(1, opU.Ntheta);
        [rho, logrho, lseInfo] = src.stable_logsumexp_rows(logn, weights, 690);
        info.logrhoMax = max(logrho(:));
        info.logrhoMin = min(logrho(:));
        info.numCappedRhoHigh = lseInfo.numCappedHigh;
    otherwise
        error('reconstruct_rho_dual_1d:UnknownMode', '未知 rho 重构方式：%s', mode);
end

rho = real(rho(:));
nU = real(nU);
end
