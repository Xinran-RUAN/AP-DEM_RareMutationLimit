function [scale, info] = mass_balance_correction_1d(rhoTrial, Kx, op, epsVal, dt, prevMass, prevRhs, par)
%MASS_BALANCE_CORRECTION_1D 用总质量平衡关系计算一步质量修正比例。
%
% 设 provisional rho 为 rhoTilde，修正后 rho^{n+1}=lambda*rhoTilde。
% 定义：
%   Mtilde = int rhoTilde dx,
%   A      = int K*rhoTilde dx,
%   B      = int rhoTilde^2 dx.
%
% backward-euler:
%   eps*(lambda*Mtilde - M^n)/dt = lambda*A - lambda^2*B.
%
% trapezoid:
%   eps*(lambda*Mtilde - M^n)/dt = 0.5*(R^n + lambda*A - lambda^2*B).
%
% 返回正根 lambda。若质量修正关闭或计算失败，lambda=1。

if nargin < 8 || isempty(par)
    par = struct();
end

info = struct();
info.applied = false;
info.enabled = local_mass_correction_enabled(par);
info.formula = local_mass_correction_formula(par);
info.scale = 1;
info.reason = 'disabled';
info.prevMass = prevMass;
info.prevReactionIntegral = prevRhs;
info.dt = dt;
info.eps = epsVal;

stats = src.mass_balance_stats_1d(rhoTrial, Kx, op);
info.trialMass = stats.mass;
info.trialReactionIntegral = stats.reactionIntegral;
info.trialA = stats.A;
info.trialB = stats.B;
info.trialRhoMin = stats.rhoMin;
info.trialRhoMax = stats.rhoMax;

scale = 1;

if ~info.enabled
    return;
end
if ~(isfinite(dt) && dt > 0 && isfinite(epsVal) && epsVal > 0 && ...
        isfinite(prevMass) && prevMass >= 0 && isfinite(prevRhs))
    info.reason = 'invalid previous mass/rhs or dt/eps';
    return;
end
if ~(isfinite(stats.mass) && stats.mass >= 0 && isfinite(stats.A) && isfinite(stats.B) && stats.B >= 0)
    info.reason = 'invalid trial stats';
    return;
end

switch lower(info.formula)
    case {'be','backward-euler','backward_euler','euler','implicit','balance-be','balance-backward-euler'}
        qa = dt * stats.B;
        qb = epsVal * stats.mass - dt * stats.A;
        qc = -epsVal * prevMass;
        info.formula = 'backward-euler';
    case {'trap','trapezoid','trapezoidal','cn','crank-nicolson','balance-trapezoid','balance-trap'}
        qa = 0.5 * dt * stats.B;
        qb = epsVal * stats.mass - 0.5 * dt * stats.A;
        qc = -(epsVal * prevMass + 0.5 * dt * prevRhs);
        info.formula = 'trapezoid';
    otherwise
        error('mass_balance_correction_1d:UnknownFormula', '未知 massCorrectionFormula: %s', info.formula);
end

info.quadraticA = qa;
info.quadraticB = qb;
info.quadraticC = qc;

scale = local_positive_root(qa, qb, qc);

minScale = local_get_num(par, 'massCorrectionMinScale', 0);
maxScale = local_get_num(par, 'massCorrectionMaxScale', Inf);
if isfinite(minScale) && minScale > 0 && scale < minScale
    scale = minScale;
    info.reason = 'clipped to minScale';
elseif isfinite(maxScale) && maxScale > 0 && scale > maxScale
    scale = maxScale;
    info.reason = 'clipped to maxScale';
else
    info.reason = 'ok';
end

if ~(isfinite(scale) && scale > 0)
    scale = 1;
    info.reason = 'invalid scale fallback to 1';
    return;
end

info.applied = true;
info.scale = scale;
info.correctedMass = scale * stats.mass;
info.correctedA = scale * stats.A;
info.correctedB = scale^2 * stats.B;
info.correctedReactionIntegral = scale * stats.A - scale^2 * stats.B;
end

function tf = local_mass_correction_enabled(par)
% 默认字段名：massCorrectionDuringSolve。也兼容 massCorrection='none'/'off'。
tf = false;
if isstruct(par) && isfield(par, 'massCorrectionDuringSolve') && ~isempty(par.massCorrectionDuringSolve)
    tf = logical(par.massCorrectionDuringSolve);
end
if isstruct(par) && isfield(par, 'massCorrection') && ~isempty(par.massCorrection)
    mode = lower(char(par.massCorrection));
    if any(strcmp(mode, {'none','off','false','0','no'}))
        tf = false;
    elseif any(strcmp(mode, {'on','true','1','yes','balance','balance-be','balance-trapezoid','trapezoid','backward-euler'}))
        tf = true;
    end
end
end

function formula = local_mass_correction_formula(par)
formula = 'trapezoid';
if isstruct(par) && isfield(par, 'massCorrectionFormula') && ~isempty(par.massCorrectionFormula)
    formula = char(par.massCorrectionFormula);
elseif isstruct(par) && isfield(par, 'massCorrection') && ~isempty(par.massCorrection)
    mode = lower(char(par.massCorrection));
    if any(strcmp(mode, {'balance-be','be','backward-euler'}))
        formula = 'backward-euler';
    elseif any(strcmp(mode, {'balance-trapezoid','trapezoid','trap'}))
        formula = 'trapezoid';
    end
end
end

function x = local_positive_root(a, b, c)
if abs(a) <= realmin
    if abs(b) <= realmin
        x = NaN;
    else
        x = -c / b;
    end
    return;
end
D = b*b - 4*a*c;
if D < 0 && D > -100*eps(max(abs([b*b,4*a*c,1])))
    D = 0;
end
if D < 0 || ~isfinite(D)
    x = NaN;
    return;
end
sD = sqrt(D);
% 二次方程通常 c<=0，因此正根唯一。根据 b 的符号选更稳定表达式。
if b >= 0
    denom = sD + b;
    if abs(denom) <= realmin
        x = (-b + sD) / (2*a);
    else
        x = -2*c / denom;
    end
else
    x = (-b + sD) / (2*a);
end
end

function v = local_get_num(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
end
