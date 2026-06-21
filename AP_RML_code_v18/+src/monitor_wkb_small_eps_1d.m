function info = monitor_wkb_small_eps_1d(W, u, Wold, uold, eps, dt, op, par, ampInfo, dtInfo, steadyInfo)
%MONITOR_WKB_SMALL_EPS_1D 小 epsilon 病态来源监测。
%   这个监测函数的目的不是改变算法，而是告诉我们：
%   - 是否相邻相位差除以 eps 过大；
%   - 峰宽是否已经小于 theta 网格；
%   - 积分因子 q 是否过大；
%   - W 是否出现负值/非有限值；
%   - residual 不下降更像是相位问题、振幅问题还是密度峰欠解析。
if nargin < 9 || isempty(ampInfo), ampInfo = struct(); end
if nargin < 10 || isempty(dtInfo), dtInfo = struct(); end
if nargin < 11 || isempty(steadyInfo), steadyInfo = struct(); end

info = struct();
up = u(:).';
duNeighbor = diff([up, up(1)]);
info.maxNeighborJumpOverEps = max(abs(duNeighbor)) / max(eps, realmin);
info.maxStepGaugeDiffOverEps = max(abs((u(:)-max(u(:))) - (uold(:)-max(uold(:))))) / max(eps, realmin);

% 峰宽估计：width ≈ sqrt(eps/(-u''(theta_*)))。
[~, kstar] = max(up);
ddu = op.Ltheta * u(:);
curv = -real(ddu(kstar));
if curv > 0
    width = sqrt(eps / curv);
    info.peakWidth = width;
    info.peakWidthOverGrid = width / op.dtheta;
else
    info.peakWidth = Inf;
    info.peakWidthOverGrid = Inf;
end

info.WMin = min(real(W(:)));
info.WMax = max(real(W(:)));
info.numNegativeW = nnz(real(W(:)) < 0);
info.numNonfiniteW = nnz(~isfinite(W(:)));
info.dtOverEps = dt / max(eps, realmin);

% 搬运振幅更新中的关键病态指标。
info.logQMax = get_field(ampInfo, 'logQMax', NaN);
info.logQMin = get_field(ampInfo, 'logQMin', NaN);
info.ifGaugeShift = get_field(ampInfo, 'ifGaugeShift', NaN);
info.numRhsLogClippedHigh = get_field(ampInfo, 'numRhsLogClippedHigh', NaN);
info.numRhsLogClippedLow = get_field(ampInfo, 'numRhsLogClippedLow', NaN);
info.matrixCondEst = get_field(ampInfo, 'matrixCondEst', NaN);
info.matrixRcondEst = get_field(ampInfo, 'matrixRcondEst', NaN);
info.phaseCFL = get_field(dtInfo, 'phaseCFL', NaN);
info.res_phase_steady = get_field(steadyInfo, 'res_phase_steady', NaN);
info.res_density_steady = get_field(steadyInfo, 'res_density_steady', NaN);
info.res_amp_w_steady = get_field(steadyInfo, 'res_amp_w_steady', NaN);

% 给出一个粗略中文提示，写入 history；进度行中也可以显示部分量。
th = local_get_num(par, 'monitorWarnThreshold', 50);
notes = {};
if isfinite(info.maxNeighborJumpOverEps) && info.maxNeighborJumpOverEps > th
    notes{end+1} = '相邻u差/eps很大'; %#ok<AGROW>
end
if isfinite(info.logQMax) && max(abs([info.logQMax, info.logQMin])) > th
    notes{end+1} = '相位积分因子q很大'; %#ok<AGROW>
end
if isfinite(info.peakWidthOverGrid) && info.peakWidthOverGrid < 1.5
    notes{end+1} = '密度峰宽小于约1.5个theta网格'; %#ok<AGROW>
end
if info.numNegativeW > 0
    notes{end+1} = 'W出现负值'; %#ok<AGROW>
end
if isempty(notes)
    info.noteCode = 'OK';
else
    info.noteCode = strjoin(notes, '; ');
end
end

function val = get_field(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end

function val = local_get_num(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end
