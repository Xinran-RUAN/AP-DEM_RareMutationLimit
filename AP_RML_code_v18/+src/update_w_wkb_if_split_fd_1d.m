function [Wnew, uWork, info] = update_w_wkb_if_split_fd_1d(W, uOld, uRaw, D, Kx, rho, H, eps, dt, op, par)
%UPDATE_W_WKB_IF_SPLIT_FD_1D  WKB 框架内的小 epsilon 稳健振幅更新。
%
%  本函数始终推进 WKB 变量 (W,u)，没有回到直接求解 n 方程。
%
%  小 epsilon 下上一版爆炸的主要来源有三个：
%    1) u 的常数 gauge 漂移，使 Laplace 重构中的 exp(u*/eps) 上溢；
%    2) full IF 或不平衡 IF 使 q=exp(logq) 的整体尺度过大；
%    3) dt/eps 很大时，反应快层导致 W 振幅出现极端缩放。
%
%  本版本的处理：
%    A. phase-peak gauge：把 uRaw 的二次插值峰值规范为 0，避免 u*/eps 上溢。
%       这个常数平移通过 q 的反向缩放补偿，不改变重构密度。
%    B. H-only IF：默认只抵消相位方程中的 -H 与振幅方程中的 +H。
%    C. exact logistic row scaling：对每个 x_j 的非局部竞争快层作解析缩放，
%       然后仍然求解 W 的线性隐式扩散/输运步。
%    D. 只要检测到 RHS log 过大、线性系统非有限、或 W 超界，就设置
%       info.needsRetry=true，让 run_wkb_1d 缩小本步 dt 后重试。

if nargin < 11 || isempty(par)
    par = struct();
end
Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;
info = struct();
info.needsRetry = false;
info.retryReason = '';

reactionMode = lower(local_get_string(par, 'reactionDiscretization', 'logistic-exact'));
ifMode = lower(local_get_string(par, 'ifMode', 'h-only'));

%-------------------------
% 1. 积分因子 logQ0。
%-------------------------
switch ifMode
    case {'h-only','h','hamiltonian-only','eigenvalue-only'}
        logQ0 = dt * real(H(:).') / eps;
        info.ifModeUsed = 'h-only';
    case {'full','full-phase'}
        logQ0 = (real(uOld(:).') - real(uRaw(:).')) / eps;
        info.ifModeUsed = 'full-phase';
    case {'none','off'}
        logQ0 = zeros(1,Ktheta);
        info.ifModeUsed = 'none';
    otherwise
        error('Unknown ifMode "%s".', ifMode);
end
info.qLog0Max = max(logQ0(:));
info.qLog0Min = min(logQ0(:));
info.qLog0Range = info.qLog0Max - info.qLog0Min;

%-------------------------
% 2. 选择 u^{n+1} 的常数 gauge。
%    默认 phase-peak：二次插值峰值为 0。这个选择比“让 Wseed 最大值为 10”
%    更适合小 epsilon，因为它从源头避免 rho 重构里的 exp(u*/eps) 上溢。
%-------------------------
gaugeStrategy = lower(local_get_string(par, 'ifGaugeStrategy', 'phase-peak'));
peakRaw = src.phase_peak_info_1d(uRaw, op);
switch gaugeStrategy
    case {'phase-peak','peak','max-phase','max-u'}
        cGauge = -peakRaw.uStar;
        gaugeReason = 'phase-peak';
    case {'node-max','nodal-max'}
        cGauge = -max(real(uRaw(:)));
        gaugeReason = 'node-max';
    case {'balanced','rhs-balance'}
        WposTmp = max(real(W), realmin);
        logSeed0Tmp = log(WposTmp) + repmat(logQ0, nx, 1);
        target = local_get_num(par, 'ifRhsLogTarget', log(local_get_num(par, 'ifSeedTargetValue', 10)));
        cGauge = eps * (max(logSeed0Tmp(:)) - target);
        gaugeReason = 'rhs-balance';
    case {'none','off'}
        cGauge = 0;
        gaugeReason = 'none';
    otherwise
        error('Unknown ifGaugeStrategy "%s".', gaugeStrategy);
end
uWork = real(uRaw) + cGauge;
peakWork = src.phase_peak_info_1d(uWork, op);

% q 需要补偿 u 的常数平移：q*exp(uWork/eps)=q0*exp(uRaw/eps)。
logQ = logQ0 - cGauge/eps;
info.ifGaugeShift = cGauge;
info.ifGaugeStrategyUsed = gaugeReason;
info.uPeakBeforeGauge = peakRaw.uStar;
info.uPeakAfterGauge = peakWork.uStar;
info.uNodeMaxAfterGauge = max(uWork(:));
info.qLogMax = max(logQ(:));
info.qLogMin = min(logQ(:));
info.logQMax = info.qLogMax;
info.logQMin = info.qLogMin;
info.qLogAbsMax = max(abs([info.qLogMax, info.qLogMin]));
info.qLogRange = info.qLogMax - info.qLogMin;

%-------------------------
% 3. 构造 Wseed=q W。只在 log 域中先检查尺度。
%-------------------------
Wpos = max(real(W), realmin);
logW = log(Wpos);
logSeed = logW + repmat(logQ, nx, 1);
info.seedLogMax = max(logSeed(:));
info.seedLogMin = min(logSeed(:));
info.seedLogRange = info.seedLogMax - info.seedLogMin;

% 这些阈值不是数学稳定性条件，而是双精度 W 变量可表示性的保护。
% v8 中把 seedLogRange 稍大就作为拒步条件，导致小 epsilon 下经常卡死。
% v10 默认只在真正接近双精度溢出时才拒步；较大的 theta 振幅跨度先记录
% 为诊断量，不再自动认为格式失败。
seedAbsMax = local_get_num(par, 'ifSeedLogAbsMax', 600);
seedRangeMax = local_get_num(par, 'ifSeedLogRangeMax', 500);
rejectOnSeed = local_get_bool(par, 'rejectOnSeedLogRange', false);
if rejectOnSeed && (info.seedLogMax > seedAbsMax || info.seedLogMin < -seedAbsMax || info.seedLogRange > seedRangeMax)
    info.needsRetry = true;
    info.retryReason = sprintf('WKB-IF RHS log 尺度过大: max=%.3g, min=%.3g, range=%.3g', ...
        info.seedLogMax, info.seedLogMin, info.seedLogRange);
end

% 为了避免 exp 溢出，这里仅在非常极端时截断。若发生截断，history 中会记录。
logClip = local_get_num(par, 'ifRhsLogHardClip', 650);
logSeedClip = min(max(logSeed, -logClip), logClip);
Wseed = exp(logSeedClip);
info.numRhsLogClippedHigh = nnz(logSeed > logClip);
info.numRhsLogClippedLow = nnz(logSeed < -logClip);
info.qClipCount = info.numRhsLogClippedHigh + info.numRhsLogClippedLow;
info.qClipFraction = info.qClipCount / max(numel(logSeed),1);

%-------------------------
% 4. exact logistic row scaling，仍然作用在 W 上。
%-------------------------
useExactReaction = any(strcmpi(reactionMode, {'logistic-exact','exact-logistic','exact','ap-reaction','patankar'}));
if strcmpi(reactionMode, 'patankar')
    info.reactionModeRequested = 'patankar';
    reactionMode = 'logistic-exact';
end
if useExactReaction
    [scaleX, rhoAfterReact, reactInfo] = local_logistic_reaction_scale(rho(:), Kx(:), eps, dt, par);
    Wseed = bsxfun(@times, Wseed, scaleX(:));
    info.rhoReactionBeforeMax = max(rho(:));
    info.rhoReactionAfterMax = max(rhoAfterReact(:));
    info.reactionScaleMax = max(scaleX(:));
    info.reactionScaleMin = min(scaleX(:));
    info.reactionModeUsed = reactionMode;
    info.reactionScaleClipped = reactInfo.numClipped;
else
    info.reactionModeUsed = reactionMode;
    info.rhoReactionBeforeMax = max(rho(:));
    info.rhoReactionAfterMax = NaN;
    info.reactionScaleMax = 1;
    info.reactionScaleMin = 1;
    info.reactionScaleClipped = 0;
end

% 若 logistic 缩放也过大，要求重试。通常这说明 rho 重构极小或 dt 太大。
scaleLimit = local_get_num(par, 'reactionScaleRetryLimit', 1e4);
if info.reactionScaleMax > scaleLimit || info.reactionScaleMin < 1/scaleLimit
    info.needsRetry = true;
    if isempty(info.retryReason)
        info.retryReason = sprintf('reaction scale 过大: [%.3g, %.3g]', info.reactionScaleMin, info.reactionScaleMax);
    end
end

%-------------------------
% 5. 剩余 W 方程的线性隐式步。
%-------------------------
ddu = (op.Ltheta * uWork(:)).';
rCurvMat = repmat((-2*eps*ddu(:).'), nx, 1);
rCurv = rCurvMat(:);
rPlus = max(rCurv, 0);
rMinus = max(-rCurv, 0);

b = eps * Wseed(:);
DblockLx = kron(spdiags(D(:),0,Ktheta,Ktheta), op.Lx);
Lth = op.Ltheta_big;
Tmat = src.upwind_transport_matrix(uWork, op);
A = eps*speye(N) - dt*DblockLx - dt*eps^2*Lth;
if ~isfield(par,'transportImplicit') || par.transportImplicit
    A = A + dt*2*eps*Tmat;
else
    b = b - dt*2*eps*(Tmat*Wseed(:));
end
A = A + dt*spdiags(rMinus,0,N,N);
b = b + dt*(rPlus(:).*Wseed(:));

info.matrixDiagMin = min(abs(diag(A)));
info.matrixDiagMax = max(abs(diag(A)));
info.matrixCondEst = NaN;
info.matrixRcondEst = NaN;
if local_get_bool(par, 'monitorMatrixCondition', false)
    try
        info.matrixCondEst = condest(A);
        info.matrixRcondEst = 1 / info.matrixCondEst;
    catch
        info.matrixCondEst = NaN;
        info.matrixRcondEst = NaN;
    end
end

Wvec = A \ b;
Wnew = reshape(real(Wvec), nx, Ktheta);
info.numNonfiniteW = nnz(~isfinite(Wnew));
info.minWBeforeClip = min(Wnew(:));
info.maxWBeforeClip = max(Wnew(:));

if info.numNonfiniteW > 0
    info.needsRetry = true;
    if isempty(info.retryReason), info.retryReason = 'W 线性求解产生非有限值'; end
    Wnew(~isfinite(Wnew)) = 0;
end

if local_get_bool(par, 'clipNegativeW', true)
    floorW = local_get_num(par, 'WFloor', 0);
    info.numNegativeW = nnz(Wnew < floorW);
    Wnew = max(Wnew, floorW);
else
    info.numNegativeW = nnz(Wnew < 0);
end

WmaxAllowed = local_get_num(par, 'WRetryMax', 1e80);
if max(abs(Wnew(:))) > WmaxAllowed
    info.needsRetry = true;
    if isempty(info.retryReason), info.retryReason = 'W 超过 WRetryMax'; end
end

info.rIfMax = max(rCurv);
info.rIfMin = min(rCurv);
info.rIfPlusMax = max(max(rCurv),0);
info.rIfAbsMax = max(abs(rCurv));
info.WmaxAfter = max(Wnew(:));
info.WminAfter = min(Wnew(:));
end

function [scale, rhoNew, info] = local_logistic_reaction_scale(rho, K, eps, dt, par)
% 对 eps*rho_t=(K-rho)rho 的逐点精确解给出缩放 rhoNew/rho。
info = struct('numClipped',0);
rho = max(real(rho(:)), 0);
K = max(real(K(:)), 0);
rhoFloor = local_get_num(par, 'rhoReactionFloor', 1e-300);
rhoSafe = max(rho, rhoFloor);
rhoNew = zeros(size(rhoSafe));

for j = 1:numel(rhoSafe)
    Kj = K(j);
    rj = rhoSafe(j);
    if Kj > 0
        a = Kj * dt / max(eps, realmin);
        if a > 80
            rhoNew(j) = Kj;  % 快层已经充分松弛到 carrying capacity
        else
            ea = exp(-a);
            rhoNew(j) = Kj / (1 + (Kj/rj - 1)*ea);
        end
    else
        rhoNew(j) = rj / (1 + dt*rj/max(eps,realmin));
    end
end
scale = rhoNew ./ rhoSafe;

sMax = local_get_num(par, 'reactionScaleMaxAllowed', 1e8);
sMin = local_get_num(par, 'reactionScaleMinAllowed', 1e-8);
n1 = 0; n2 = 0;
if isfinite(sMax) && sMax > 0
    n1 = nnz(scale > sMax);
    scale = min(scale, sMax);
end
if isfinite(sMin) && sMin > 0
    n2 = nnz(scale < sMin);
    scale = max(scale, sMin);
end
info.numClipped = n1 + n2;
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end

function tf = local_get_bool(s, name, defaultValue)
tf = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end

function val = local_get_num(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end
