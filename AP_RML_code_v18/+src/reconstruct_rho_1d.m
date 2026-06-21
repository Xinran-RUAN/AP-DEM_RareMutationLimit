function [rho, n, info] = reconstruct_rho_1d(W, u, eps, op, mode)
%RECONSTRUCT_RHO_1D  由 WKB 变量重构 n 和 rho。
%
%  n(x,theta) = W(x,theta) * exp(u(theta)/eps),
%  rho(x)     = int n(x,theta) dtheta.
%
%  注意：本函数只是 WKB 变量的重构步骤，不是原始密度方程的求解器。
%  小 epsilon 下 direct-log 主要用于诊断；正式 rho 推荐 laplace-hybrid。

if nargin < 5 || isempty(mode)
    mode = 'direct-log';
end
info = struct('mode', mode, 'usedFallback', false, 'thetaStar', NaN, 'dduStar', NaN, ...
              'numNonpositiveW', 0, 'lognMax', NaN, 'logrhoMax', NaN, ...
              'numCappedNHigh', 0, 'numCappedRhoHigh', 0);

switch lower(mode)
    case {'direct','direct-log'}
        Wreal = real(W);
        bad = (~isfinite(Wreal)) | (Wreal <= 0);
        info.numNonpositiveW = nnz(bad);
        Wpos = max(Wreal, realmin);
        logW = log(Wpos);
        logn = bsxfun(@plus, logW, real(u(:).')/eps);
        info.lognMax = max(logn(:));
        info.lognMin = min(logn(:));
        info.numNonfiniteLogN = nnz(~isfinite(logn(:)));

        % n 用于图像和诊断，不用于推进。若 logn 过大，截断并记录。
        logNCap = 690;
        info.numCappedNHigh = nnz(logn > logNCap);
        n = exp(min(logn, logNCap));
        n(~isfinite(n)) = 0;

        weights = op.dtheta * ones(1, op.Ntheta);
        [rho, logrho, lseInfo] = src.stable_logsumexp_rows(logn, weights, 690);
        info.logrhoMax = max(logrho(:));
        info.logrhoMin = min(logrho(:));
        info.numCappedRhoHigh = lseInfo.numCappedHigh;

    case 'laplace'
        [rho, n, info] = src.reconstruct_laplace_1d(W, u, eps, op, false);
    case {'laplace-hybrid','hybrid-laplace'}
        [rho, n, info] = src.reconstruct_laplace_1d(W, u, eps, op, true);
    otherwise
        error('Unknown rho reconstruction mode "%s".', mode);
end
rho = real(rho(:));
n = real(n);
end
