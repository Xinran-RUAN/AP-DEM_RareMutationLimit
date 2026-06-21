function info = compute_wkb_steady_residual_1d(W, u, D, Kx, rho, H, eps, op, par, variant)
%COMPUTE_WKB_STEADY_RESIDUAL_1D 计算 WKB 稳态方程残差。
%
%  这个残差直接检查 WKB 稳态系统，而不是把主求解器改为 n 方程，也不是
%  简单检查单步差分。它帮助判断 residual 不下降究竟来自相位方程还是 W 振幅方程。
%  其中 density residual 仅作为 gauge-invariant 诊断量保存，不参与主推进。

if nargin < 10 || isempty(variant)
    variant = local_get_string(par, 'amplitudeVariant', 'split-if');
end

Hhat = src.numerical_hamiltonian(u, op.dtheta, par.phaseHamiltonian, par.lfAlpha);
ddu = (op.Ltheta * u(:)).';
phaseVec = Hhat(:).' - eps*ddu(:).' + H(:).';
phaseDen = 1 + max(abs(H(:))) + max(abs(Hhat(:)));
resPhase = max(abs(phaseVec)) / max(phaseDen, realmin);

Ktheta = op.Ntheta;
DblockLx = kron(spdiags(D(:),0,Ktheta,Ktheta), op.Lx);
Lth = op.Ltheta_big;
Tmat = src.upwind_transport_matrix(u, op);
Wvec = W(:);

% 原始 split-WKB 稳态振幅方程：
% -D Lx W - eps^2 Ltheta W + 2eps DthetaF(W,u)
% = W*(K-rho+H-2eps*u_thetatheta)。
extra = -2*eps*ddu;
rSplit = src.reaction_vector_1d(Kx, rho, H, extra, op);
ampVec = -DblockLx*Wvec - eps^2*(Lth*Wvec) + 2*eps*(Tmat*Wvec) - rSplit(:).*Wvec;
ampScale = 1 + max(abs(rSplit(:).*Wvec)) + max(abs(DblockLx*Wvec)) + ...
           max(abs(eps^2*(Lth*Wvec))) + max(abs(2*eps*(Tmat*Wvec)));
resAmp = max(abs(ampVec)) / max(ampScale, realmin);

% 仅诊断：重构密度层面的稳态残差。小 eps 粗 theta 网格下它会反映峰形欠解析，
% 不代表 WKB phase 的 fittest trait 定位一定失败。
[~, n] = src.reconstruct_rho_1d(W, u, eps, op, 'direct-log');
nvec = n(:);
rBase = Kx(:) - rho(:);
rBaseVec = repmat(rBase, Ktheta, 1);
densVec = -DblockLx*nvec - eps^2*(Lth*nvec) - rBaseVec(:).*nvec;
densScale = 1 + max(abs(rBaseVec(:).*nvec)) + max(abs(DblockLx*nvec)) + max(abs(eps^2*(Lth*nvec)));
resDens = max(abs(densVec)) / max(densScale, realmin);

% 若出现 Inf/NaN，不能把 residual 误判为 0；直接返回 Inf 并让主程序报警。
if ~isfinite(resPhase), resPhase = Inf; end
if ~isfinite(resAmp),   resAmp   = Inf; end
if ~isfinite(resDens),  resDens  = Inf; end

info = struct();
info.res_phase_steady = real(resPhase);
info.res_amplitude_steady = real(resAmp);
info.res_amp_w_steady = real(resAmp);
info.res_density_steady = real(resDens);
info.res_wkb_steady = max(info.res_phase_steady, info.res_amplitude_steady);
info.phase_raw_max = max(abs(phaseVec));
info.amplitude_raw_max = max(abs(ampVec));
info.density_raw_max = max(abs(densVec));
if ~isfinite(info.phase_raw_max), info.phase_raw_max = Inf; end
if ~isfinite(info.amplitude_raw_max), info.amplitude_raw_max = Inf; end
if ~isfinite(info.density_raw_max), info.density_raw_max = Inf; end
info.variantChecked = variant;
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end
