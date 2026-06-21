function [dtStep, info] = select_wkb_time_step_1d(par, t, W, u, D, Kx, rho, H, op, variant)
%SELECT_WKB_TIME_STEP_1D 可选自适应时间步。
%   split-if 默认限制相位 Hamilton-Jacobi 显式 CFL；对 WKB-IF 还可限制
%   IF 指数在 theta 方向的相对跨度 max(log q)-min(log q)，避免 W 振幅
%   出现不可表示的指数分离。这不是解 n 方程，只是 WKB 变量的数值保护。

if nargin < 10 || isempty(variant)
    variant = 'split';
end
baseDt = get_num(par, 'dt', 1e-3);
T = get_num(par, 'T', Inf);
eps = get_num(par, 'eps', 1);
dtMax = baseDt;
if isfield(par, 'dtMax') && ~isempty(par.dtMax) && isnumeric(par.dtMax) && isscalar(par.dtMax)
    if isfinite(double(par.dtMax)) && double(par.dtMax) > 0
        dtMax = min(baseDt, double(par.dtMax));
    end
end
remaining = T - t;
if isfinite(remaining)
    dtMax = min(dtMax, max(remaining, 0));
end

useAdaptive = false;
if isfield(par, 'adaptiveTimeStep') && ~isempty(par.adaptiveTimeStep)
    useAdaptive = logical(par.adaptiveTimeStep);
elseif isfield(par, 'useAdaptiveDt') && ~isempty(par.useAdaptiveDt)
    useAdaptive = logical(par.useAdaptiveDt);
end

info = src.time_step_stability_info_1d(W, u, Kx, rho, H, eps, max(dtMax,realmin), op, par, variant);
info.dtSuggested = dtMax;
info.dtLimited = false;

if useAdaptive
    strategy = lower(get_string(par, 'adaptiveStrategy', 'phase-cfl'));
    if strcmp(strategy, 'auto')
        if any(strcmp(lower(char(variant)), {'split-if','wkb-if','wkb-if-split'}))
            strategy = 'phase-cfl';
        else
            strategy = 'combined';
        end
    end
    dtCand = dtMax;
    switch strategy
        case {'phase-cfl','cfl','hj-cfl'}
            dtCand = min(dtCand, phase_cfl_dt(par, info, op));
            dtCand = min(dtCand, if_log_range_dt(par, eps, H, variant));
        case {'eps-source','source','old'}
            dtCand = min(dtCand, eps_source_dt(par, eps, info));
        case {'combined','both'}
            dtCand = min(dtCand, phase_cfl_dt(par, info, op));
            dtCand = min(dtCand, if_log_range_dt(par, eps, H, variant));
            dtCand = min(dtCand, eps_source_dt(par, eps, info));
        case {'none','off'}
            dtCand = dtMax;
        otherwise
            error('Unknown adaptiveStrategy "%s".', strategy);
    end
    info.adaptiveStrategyUsed = strategy;
    info.dtSuggested = dtCand;
    if isfinite(dtCand) && dtCand > 0
        info.dtLimited = dtCand < dtMax;
        dtMax = min(dtMax, dtCand);
    end
end

dtMin = get_num(par, 'dtMin', 0);
if isfinite(dtMin) && dtMin > 0 && dtMax < dtMin && remaining > dtMin
    dtMax = dtMin;
end

dtStep = max(dtMax, 0);
info.dt = dtStep;
info.dtOverEps = dtStep / eps;
info.lambdaR = dtStep / eps * info.rplusMax;
info.lambdaH = dtStep / eps * info.HabsMax;
info.phaseCFL = info.phaseSpeedMax * dtStep / op.dtheta;
info.ifLogRangeEstimate = dtStep / eps * max(max(H(:))-min(H(:)),0);
end

function dtc = phase_cfl_dt(par, info, op)
safety = get_num(par, 'phaseCflSafety', 0.45);
if ~(isfinite(safety) && safety > 0), safety = 0.45; end
if info.phaseSpeedMax <= 0 || ~isfinite(info.phaseSpeedMax)
    dtc = Inf;
else
    dtc = safety * op.dtheta / info.phaseSpeedMax;
end
end


function dtc = if_log_range_dt(par, eps, H, variant)
% WKB-IF 中即使使用 gauge 去掉整体 q，theta 间的 q 相对跨度
% exp((max H-min H)*dt/eps) 仍会直接进入 W 的 theta 方向振幅。
% 该跨度过大时，W 变量本身会指数分离，双精度线性代数无法可靠表示。
useIt = get_bool(par, 'adaptiveIfLogRange', true);
if ~useIt || ~any(strcmp(lower(char(variant)), {'wkb-if-split','wkb-if','split-if'}))
    dtc = Inf;
    return;
end
L = get_num(par, 'ifLogRangeMax', 80);
if ~(isfinite(L) && L > 0), L = 80; end
hr = max(H(:)) - min(H(:));
if ~(isfinite(hr) && hr > 0)
    dtc = Inf;
else
    dtc = L * eps / hr;
end
end

function tf = get_bool(s, name, defaultValue)
tf = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end

function dtc = eps_source_dt(par, eps, info)
eta = get_num(par, 'adaptiveSafety', 0.5);
if ~(isfinite(eta) && eta > 0), eta = 0.5; end
denom = max([1, info.rplusMax, info.HabsMax]);
dtc = eta * eps / denom;
end

function val = get_num(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end

function val = get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end
