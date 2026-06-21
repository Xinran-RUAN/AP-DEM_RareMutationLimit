function [Wout, uout, info] = stabilize_wkb_gauge_1d(W, u, eps, par)
%STABILIZE_WKB_GAUGE_1D 用 gauge 自由度控制 W 的数值量级。
%
%  该函数不改变 n = W exp(u/eps)，只改变 W 与 u 之间的常数分配。
%  在小 eps 下，这一步可以避免 W 因 max-gauge 被指数放大，从而改善 residual
%  和线性代数稳定性。

if nargin < 4 || isempty(par)
    par = struct();
end
mode = local_get_string(par, 'postGaugeMode', local_get_string(par, 'gaugeMode', 'incremental'));
info = struct('postGaugeMode', mode, 'postGaugeShift', 0);

switch lower(mode)
    case {'max','max-zero'}
        [Wout, uout, c] = src.normalize_gauge(W, u, eps, 'max');
        info.postGaugeShift = c;
    case {'none','no','incremental','balanced','if-balanced','mean','mean-increment'}
        % 默认不再强行 max(u)=0。因为前面的 choose_phase_gauge_1d 已经选择了
        % 对 IF 更新更合适的 gauge。
        Wout = W;
        uout = u;
    case {'w-balance','amplitude-balance','balance-w'}
        % 将 W 的典型大小保持在 target 附近。只用正的有限 W 估计 logW。
        target = local_get_num(par, 'targetLogW', 0.0);
        Wpos = W(isfinite(W) & W > 0);
        if isempty(Wpos)
            Wout = W; uout = u;
        else
            logW = log(Wpos(:));
            center = median(logW);
            c = eps*(target - center);
            uout = u - c;
            Wout = W * exp(c/eps);
            info.postGaugeShift = c;
        end
    otherwise
        error('Unknown postGaugeMode "%s".', mode);
end
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end

function val = local_get_num(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end
