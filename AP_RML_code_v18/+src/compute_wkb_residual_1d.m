function [residual, info] = compute_wkb_residual_1d(Wnew, unew, Wold, uold, eps, dt, op, mode)
%COMPUTE_WKB_RESIDUAL_1D Gauge-aware time residual for WKB variables.
%   mode='scaled-gauge-invariant' (default for small epsilon) uses
%
%       max( ||u^{n+1}-u^n||_gauge/dt,
%            eps*||n^{n+1}-n^n||/(dt*(1+||n^{n+1}||)) ).
%
%   The factor eps is natural because the original density equation is
%       eps*n_t = ...
%   Without this factor, a pseudo-time iteration in the rare-mutation regime
%   can report artificially large residuals even when the steady equation
%   residual is moderate.  The unscaled density time residual is still stored
%   as info.res_n.
%
%   mode='gauge-invariant' keeps the older unscaled density residual.
%   mode='legacy' uses the old W-based residual, which is not gauge-invariant.
if nargin < 8 || isempty(mode)
    mode = 'scaled-gauge-invariant';
end
mode = lower(char(mode));

unewGauge = unew(:) - max(unew(:));
uoldGauge = uold(:) - max(uold(:));
resU = max(abs(unewGauge - uoldGauge)) / dt;

[~, nnew] = src.reconstruct_rho_1d(Wnew, unew, eps, op, 'direct-log');
[~, nold] = src.reconstruct_rho_1d(Wold, uold, eps, op, 'direct-log');
resN = max(abs(nnew(:)-nold(:))) / (dt*(1 + max(abs(nnew(:)))));
resNScaled = eps * resN;
resWLegacy = max(abs(Wnew(:)-Wold(:))) / (dt*(1 + max(abs(Wnew(:)))));
if ~isfinite(resU), resU = Inf; end
if ~isfinite(resN), resN = Inf; end
if ~isfinite(resNScaled), resNScaled = Inf; end
if ~isfinite(resWLegacy), resWLegacy = Inf; end

switch mode
    case {'scaled-gauge-invariant','scaled-density','eps-scaled','scaled','wkb-steady','steady','steady-wkb'}
        residual = max(resU, resNScaled);
    case {'gauge-invariant','density','density-gauge','invariant'}
        residual = max(resU, resN);
    case {'density-only','n'}
        residual = resN;
    case {'density-scaled-only','n-scaled'}
        residual = resNScaled;
    case {'legacy','legacy-update','w','old'}
        residual = max(resU, resWLegacy);
    otherwise
        error('Unknown residualMode "%s".', mode);
end

info = struct();
info.mode = mode;
info.res_u = resU;
info.res_n = resN;
info.res_n_scaled = resNScaled;
info.res_W_legacy = resWLegacy;
if ~isfinite(residual), residual = Inf; end
info.residual = residual;
end
