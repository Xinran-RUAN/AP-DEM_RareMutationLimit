function info = time_step_stability_info_1d(W, u, Kx, rho, H, eps, dt, op, par, variant)
%TIME_STEP_STABILITY_INFO_1D 时间步和小 epsilon 病态诊断量。
%   这些量用于监测，不是严格稳定性定理。split-if 中零阶 H 主要由积分因子
%   处理，因此默认不把 dt 强制限制为 O(eps)。

if nargin < 10 || isempty(variant)
    variant = 'split';
end
if nargin < 9 || isempty(par)
    par = struct();
end
variant = lower(char(variant));
Ktheta = op.Ntheta;

switch variant
    case {'split','split-wkb'}
        ddu = (op.Ltheta * u(:)).';
        extra = -2*eps*ddu;
        r = src.reaction_vector_1d(Kx, rho, H, extra, op);
    case {'split-if','wkb-if','wkb-if-split','if-split'}
        ddu = (op.Ltheta * u(:)).';
        Hhat = src.numerical_hamiltonian(u, op.dtheta, par.phaseHamiltonian, par.lfAlpha);
        rIfMat = bsxfun(@plus, Kx(:)-rho(:), -Hhat(:).' - eps*ddu(:).');
        r = rIfMat(:);
    case {'density-compatible-semidiscrete','dc-semidiscrete','density-compatible-old','dc-old'}
        ddu = (op.Ltheta * u(:)).';
        Hhat = src.numerical_hamiltonian(u, op.dtheta, par.phaseHamiltonian, par.lfAlpha);
        extra = Hhat - eps*ddu;
        r = src.reaction_vector_1d(Kx, rho, H, extra, op);
    otherwise
        rBase = Kx(:) - rho(:);
        r = repmat(rBase, Ktheta, 1);
end

[dm, dp] = src.one_sided_derivatives(u, op.dtheta, 'first-order');
phaseSpeedMax = 2*max(abs([dm(:); dp(:)]));
neighborJump = max(abs(diff([u(:); u(1)])));

info = struct();
info.rplusMax = max(max(r), 0);
info.rabsMax = max(abs(r));
info.HabsMax = max(abs(H(:)));
info.dtOverEps = dt / eps;
info.lambdaR = dt / eps * info.rplusMax;
info.lambdaH = dt / eps * info.HabsMax;
info.ifLogRangeEstimate = dt / eps * max(max(H(:))-min(H(:)),0);
info.phaseSpeedMax = phaseSpeedMax;
info.phaseCFL = phaseSpeedMax * dt / op.dtheta;
info.neighborJumpOverEps = neighborJump / max(eps, realmin);
info.maxNeighborJumpOverEps = info.neighborJumpOverEps;
if ~isempty(W)
    info.Wmin = min(W(:));
    info.Wmax = max(W(:));
else
    info.Wmin = NaN;
    info.Wmax = NaN;
end
[~, kstar] = max(u(:));
ddu2 = op.Ltheta * u(:);
curv = -real(ddu2(kstar));
if curv > 0
    info.peakWidth = sqrt(eps / curv);
    info.peakWidthOverGrid = info.peakWidth / op.dtheta;
else
    info.peakWidth = Inf;
    info.peakWidthOverGrid = Inf;
end
end
