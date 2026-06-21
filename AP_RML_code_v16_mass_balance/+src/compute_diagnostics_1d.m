function diag = compute_diagnostics_1d(W, u, D, Kx, rho, H, eps, op, par, variant)
%COMPUTE_DIAGNOSTICS_1D 一维 WKB 标准诊断量。
%
%  合并版默认保持旧版定义：theta_wkb 取相位 u 的节点最大值。
%  若 par.useSubgridPhasePeak=true，则同时使用三点二次插值给出亚网格峰值，
%  这对小 epsilon 或 coarse theta 网格诊断更平滑，但属于 v13 可选补丁。

if nargin < 10 || isempty(variant)
    variant = par.amplitudeVariant;
end

[~, n] = src.reconstruct_rho_1d(W, u, eps, op, 'direct-log');
wx = utils.trapz_weights(op.nx, op.dx);
P = (wx.' * n);
[~, idP] = max(P);

[~, idUNode] = max(real(u(:)));
peak = src.phase_peak_info_1d(u, op);
useSubgrid = isfield(par,'useSubgridPhasePeak') && ~isempty(par.useSubgridPhasePeak) && logical(par.useSubgridPhasePeak);
if useSubgrid && isfinite(peak.thetaStar)
    thetaWKB = peak.thetaStar;
    idUForH = peak.k0;
else
    thetaWKB = op.theta(idUNode);
    idUForH = idUNode;
end

diag.theta_direct = op.theta(idP);
diag.theta_wkb_node = op.theta(idUNode);
diag.theta_wkb_subgrid = peak.thetaStar;
diag.theta_wkb = thetaWKB;
diag.theta_m = par.theta_m;
diag.err_direct_to_m = utils.periodic_distance(diag.theta_direct, par.theta_m);
diag.err_wkb_to_m = utils.periodic_distance(diag.theta_wkb, par.theta_m);
dist = utils.periodic_distance(op.theta, diag.theta_wkb);
diag.sigma_theta = sqrt(sum((dist.^2).*P) / max(sum(P), realmin));
diag.Hmin = min(H);
diag.H_at_wkb = H(idUForH);
diag.rho_min = min(rho);
diag.rho_max = max(rho);
diag.W_min = min(W(:));
diag.W_max = max(W(:));
diag.mass_total = op.dtheta * sum(P);
diag.gammaR = src.compute_remainder_gamma_1d(W, u, D, eps, op, par, variant);

diag.uPeakStar = peak.uStar;
diag.uPeakNode = peak.uNodeMax;
diag.uPeakCurvature = peak.dduStar;
diag.uPeakInterpolationCorrection = peak.uStar - peak.uNodeMax;
diag.useSubgridPhasePeak = double(useSubgrid);

if isfield(par, 'dt') && ~isempty(par.dt)
    ts = src.time_step_stability_info_1d(W, u, Kx, rho, H, eps, par.dt, op, par, variant);
    diag.lambdaR = ts.lambdaR;
    diag.lambdaH = ts.lambdaH;
    diag.rplusMax = ts.rplusMax;
    diag.HabsMax = ts.HabsMax;
    diag.phaseCFL = ts.phaseCFL;
    diag.phaseSpeedMax = ts.phaseSpeedMax;
    diag.neighborJumpOverEps = ts.neighborJumpOverEps;
end
Wpos = W(isfinite(W) & W > 0);
if isempty(Wpos)
    diag.logWRange = Inf;
else
    diag.logWRange = max(log(Wpos(:))) - min(log(Wpos(:)));
end
diag.P = P;
end
