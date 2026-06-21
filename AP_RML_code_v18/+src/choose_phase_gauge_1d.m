function [uGauge, c, info] = choose_phase_gauge_1d(uOld, uNewRaw, eps, par)
%CHOOSE_PHASE_GAUGE_1D 为新相位选择数值 gauge。
%
%  WKB 分解有常数 gauge 自由度：
%       u -> u-c,   W -> W exp(c/eps),
%  重构密度 n = W exp(u/eps) 不变。
%
%  小 eps 下，如果每步都强制 max(u)=0，则 W 可能被 exp(c/eps) 人为放大，
%  residual 也会被污染。这里的 incremental/balanced gauge 不改变选中特征
%  argmax u，只是选择一个更适合线性代数的常数平移。

if nargin < 4 || isempty(par)
    par = struct();
end
mode = local_get_string(par, 'gaugeMode', 'incremental');

switch lower(mode)
    case {'max','max-zero'}
        c = max(uNewRaw);
    case {'none','no'}
        c = 0;
    case {'incremental','balanced','if-balanced'}
        % 令 qLog=(uOld-(uNewRaw-c))/eps 的中位数约为 0，
        % 从而避免 q_k 整体过大或过小。
        c = median(uNewRaw(:) - uOld(:));
    case {'mean-increment','mean'}
        c = mean(uNewRaw(:) - uOld(:));
    otherwise
        error('Unknown gaugeMode "%s".', mode);
end

uGauge = uNewRaw - c;
qLog = (uOld(:).' - uGauge(:).')/eps;
info = struct();
info.gaugeMode = mode;
info.phaseGaugeShift = c;
info.qLogMinPred = min(qLog);
info.qLogMaxPred = max(qLog);
info.qLogAbsMaxPred = max(abs(qLog));
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end
