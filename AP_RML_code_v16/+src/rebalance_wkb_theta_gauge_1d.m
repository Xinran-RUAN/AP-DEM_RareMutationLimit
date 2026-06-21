function [Wb, ub, info] = rebalance_wkb_theta_gauge_1d(W, u, eps, op, par)
%REBALANCE_WKB_THETA_GAUGE_1D  WKB 变量的 theta 方向再分解。
%
%  目的：
%  小 epsilon 下，即使重构密度 n = W exp(u/eps) 很正常，振幅 W 也可能在
%  theta 方向积累很大的指数尺度。这会导致下一步 WKB-IF 右端 qW 的 log
%  范围很大，线性系统和 residual 诊断变得病态。这个现象不是 WENO 导数
%  精度造成的，而是 WKB 分解本身的非唯一性没有被充分利用。
%
%  关键思想：
%  对任意只依赖 theta 的函数 a(theta)，变换
%
%      u_new(theta)      = u(theta) + eps*a(theta) - C,
%      W_new(x,theta)    = W(x,theta)*exp(-a(theta) + C/eps)
%
%  精确保持
%
%      W_new exp(u_new/eps) = W exp(u/eps).
%
%  因此它不是回到原密度方程，而是在 WKB 框架内重新选择 phase 和 amplitude
%  的分配方式。这里选 a(theta) 为 log(W) 在 x 方向的代表值，使 W 在 theta
%  方向保持温和，而把 theta 方向的指数变化吸收到 u 中。由于加入到 u 中的
%  是 eps*a(theta)，它是 WKB 意义下的次阶相位修正，但能显著改善小 epsilon
%  的数值可表示性。
%
%  输入输出：
%    W,u  : 当前 WKB 变量。
%    Wb,ub: 再分解后的 WKB 变量，重构密度严格不变。
%    info : 记录 logW 范围、相位平移等诊断量。

info = struct();
Wreal = real(W);
u = real(u(:).');

useRephase = true;
if nargin >= 5 && isstruct(par) && isfield(par,'thetaRephase') && ~isempty(par.thetaRephase)
    useRephase = logical(par.thetaRephase);
end
if ~useRephase
    Wb = Wreal;
    ub = u;
    info.thetaRephaseApplied = 0;
    return;
end

%-----------------------------
% 1. 安全计算 log(W)。W 理论上应非负；若有极小负数，前面的格式/舍入可能
%    已产生污染。这里不直接报错，而是截断到 WFloor 并记录。
%-----------------------------
WFloor = local_get_num(par, 'thetaRephaseWFloor', 1e-300);
WFloor = max(WFloor, realmin);
numBad = nnz(~isfinite(Wreal(:)) | Wreal(:) <= 0);
Wpos = max(Wreal, WFloor);
logW = log(Wpos);
logW(~isfinite(logW)) = log(WFloor);

info.thetaRephaseApplied = 1;
info.thetaRephaseNumBadW = numBad;
info.logWRangeBefore = max(logW(:)) - min(logW(:));
info.logWMaxBefore = max(logW(:));
info.logWMinBefore = min(logW(:));

%-----------------------------
% 2. 取 x 方向的代表值 a(theta)。
%    median 比 mean 更抗局部异常值；如果需要更平滑，可改为 mean。
%-----------------------------
mode = lower(local_get_string(par, 'thetaRephaseStatistic', 'median'));
switch mode
    case {'median','x-median','robust'}
        a = median(logW, 1);
    case {'mean','x-mean'}
        a = mean(logW, 1);
    otherwise
        error('Unknown thetaRephaseStatistic "%s".', mode);
end

% 可选平滑，避免 a(theta) 中有网格噪声。默认做很弱的周期三点平滑。
smoothPasses = round(local_get_num(par, 'thetaRephaseSmoothPasses', 1));
for m = 1:max(0,smoothPasses)
    a = 0.25*circshift(a,[0,1]) + 0.50*a + 0.25*circshift(a,[0,-1]);
end

% 防止一次吸收过大的异常值。这个截断作用在 a 的相对范围上，不改变主要机制。
aMean = mean(a);
aRel = a - aMean;
clipA = local_get_num(par, 'thetaRephaseMaxAbsShift', 500);
if isfinite(clipA) && clipA > 0
    aRel = min(max(aRel, -clipA), clipA);
end
a = aMean + aRel;
info.thetaRephaseARange = max(a) - min(a);
info.thetaRephaseAMax = max(a);
info.thetaRephaseAMin = min(a);

%-----------------------------
% 3. 选择常数 C，使新的相位在节点意义下 max u = 0。
%    s = u/eps + a 是重构密度的 theta 指数包络。
%    令 C/eps=max s，可避免 Laplace/direct 重构中的 exp(u/eps) 上溢。
%-----------------------------
s = u/max(eps,realmin) + a;
C_over_eps = max(s);
ub = u + eps*a - eps*C_over_eps;

% 用 log 形式更新 W，避免中间 exp 溢出。
logWb = logW - repmat(a, size(Wreal,1), 1) + C_over_eps;
logHard = local_get_num(par, 'thetaRephaseLogHardClip', 690);
info.thetaRephaseNumClipHigh = nnz(logWb > logHard);
info.thetaRephaseNumClipLow = nnz(logWb < -logHard);
logWb = min(max(logWb, -logHard), logHard);
Wb = exp(logWb);

% 最后用二次插值峰值再做一个很小的常数 gauge 修正。这个修正确保
% Laplace 重构中的亚网格峰值也接近 0。
usePeakGauge = local_get_bool(par, 'thetaRephaseUsePeakGauge', true);
if usePeakGauge && nargin >= 4 && ~isempty(op)
    try
        pk = src.phase_peak_info_1d(ub, op);
        c = pk.uStar;
        % ub <- ub-c, Wb <- Wb*exp(c/eps)，仍保持 n 不变。
        cg = c / max(eps,realmin);
        if isfinite(cg)
            cgClip = local_get_num(par, 'thetaRephasePeakGaugeClip', 300);
            cgUse = min(max(cg, -cgClip), cgClip);
            ub = ub - eps*cgUse;
            Wb = Wb * exp(cgUse);
            info.thetaRephasePeakGauge = eps*cgUse;
            info.thetaRephasePeakGaugeClipped = double(abs(cgUse-cg) > 0);
        else
            info.thetaRephasePeakGauge = NaN;
            info.thetaRephasePeakGaugeClipped = 1;
        end
    catch
        info.thetaRephasePeakGauge = NaN;
        info.thetaRephasePeakGaugeClipped = 1;
    end
else
    info.thetaRephasePeakGauge = 0;
    info.thetaRephasePeakGaugeClipped = 0;
end

% 输出诊断。
Wpos2 = max(real(Wb), WFloor);
logW2 = log(Wpos2);
info.logWRangeAfter = max(logW2(:)) - min(logW2(:));
info.logWMaxAfter = max(logW2(:));
info.logWMinAfter = min(logW2(:));
info.uNodeMaxAfterRephase = max(ub(:));
info.uNodeMinAfterRephase = min(ub(:));
info.thetaRephaseCOverEps = C_over_eps;
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if nargin >= 1 && isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end

function tf = local_get_bool(s, name, defaultValue)
tf = defaultValue;
if nargin >= 1 && isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end

function val = local_get_num(s, name, defaultValue)
val = defaultValue;
if nargin >= 1 && isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end
