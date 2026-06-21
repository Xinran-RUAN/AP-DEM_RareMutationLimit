function [Wnew, unew, info] = wkb_one_step_frozenH_1d(Wold, uold, D, Kx, rhoCoef, Hcoef, eps, dt, op, par, variant)
%WKB_ONE_STEP_FROZENH_1D  在给定 rho,H 系数下推进一个 WKB 时间步。
%
%  这个函数是 v11 新增的公共“一步算子”。它只做一件事：
%  在当前给定的非局部密度 rhoCoef 和有效 Hamiltonian Hcoef 下，
%  从旧状态 (Wold,uold) 计算一个试探新状态 (Wnew,unew)。
%
%  注意：
%    1. 主变量始终是 W 和 u；这里没有把主问题改写成直接求 n 的方程。
%    2. rhoCoef,Hcoef 可以取旧时刻值，也可以来自全隐式 Picard 迭代中的
%       当前猜测值。因此同一个一步算子可以服务于 frozen-H、全隐式耦合和
%       projective/micro-macro 三种时间积分器。
%    3. Hcoef 必须已经完成与相位方程一致的 gauge 处理，例如 H-min(H)。
%
%  输出 info 中会记录振幅更新的 q 指数、矩阵病态、是否建议拒步等信息。

if nargin < 11 || isempty(variant)
    variant = 'wkb-if-split';
end
variant = lower(char(variant));
info = struct();
info.failed = false;
info.failureReason = '';

%-------------------------
% 1. 相位步。
%    这里仍采用原代码中的半隐式 HJ 时间离散：Hamiltonian 显式、
%    epsilon*u_{theta theta} 隐式。全隐式耦合方案中，Hcoef 由当前
%    Picard 猜测更新，因此收敛后等价于使用 H^{n+1}。
%-------------------------
uRaw = src.step_phase_1d(uold, Hcoef, eps, dt, op, par.phaseHamiltonian, par.lfAlpha);

%-------------------------
% 2. 振幅步。不同 variant 都更新 W，而不是直接求解 n。
%-------------------------
ampInfo = struct();
gaugeInfo = struct();
postGaugeInfo = struct();
try
    switch variant
        case {'wkb-if-split','wkb-if','balanced-wkb-if'}
            [Wnew, unew, ampInfo] = src.update_w_wkb_if_split_fd_1d( ...
                Wold, uold, uRaw, D, Kx, rhoCoef, Hcoef, eps, dt, op, par);
            gaugeInfo = struct('phaseGaugeShift', get_num(ampInfo,'ifGaugeShift',0));

        case {'split-if','balanced-split','if-split'}
            [uAmp, ~, gaugeInfo] = src.choose_phase_gauge_1d(uold, uRaw, eps, par);
            [Wstep, ampInfo] = src.update_w_split_if_fd_1d( ...
                Wold, uold, uAmp, D, Kx, rhoCoef, Hcoef, eps, dt, op, par);
            [Wnew, unew, postGaugeInfo] = src.stabilize_wkb_gauge_1d(Wstep, uAmp, eps, par);

        case {'split','split-wkb'}
            [Wstep, ampInfo] = src.update_w_split_fd_1d( ...
                Wold, uRaw, D, Kx, rhoCoef, Hcoef, eps, dt, op, par);
            [Wnew, unew] = src.normalize_gauge(Wstep, uRaw, eps, 'max');

        case {'density-compatible-semidiscrete','dc-semidiscrete','density-compatible-old','dc-old'}
            [Wstep, ampInfo] = src.update_w_density_compatible_fd_1d( ...
                Wold, uRaw, D, Kx, rhoCoef, Hcoef, eps, dt, op, par);
            [Wnew, unew] = src.normalize_gauge(Wstep, uRaw, eps, 'max');

        case {'density-compatible-if','dc-if','density-if','integrating-factor'}
            % 这里保留旧入口以便对照；小 epsilon 下不推荐作为默认，因为它可能
            % 形成 exp((u_{k+1}-u_k)/eps) 的指数比值矩阵。
            [Wstep, ampInfo] = src.update_w_density_if_fd_1d( ...
                Wold, uold, uRaw, D, Kx, rhoCoef, eps, dt, op, par);
            [Wnew, unew] = src.normalize_gauge(Wstep, uRaw, eps, 'max');

        otherwise
            error('Unknown amplitude variant "%s".', variant);
    end
catch ME
    Wnew = Wold;
    unew = uold;
    info.failed = true;
    info.failureReason = ['一步振幅/相位更新出错: ' ME.message];
    info.ampInfo = ampInfo;
    return;
end

%-------------------------
% 3. 基础安全检查。
%-------------------------
if isstruct(ampInfo) && isfield(ampInfo,'needsRetry') && ampInfo.needsRetry
    info.failed = true;
    info.failureReason = ampInfo.retryReason;
end
if any(~isfinite(Wnew(:))) || any(~isfinite(unew(:)))
    info.failed = true;
    if isempty(info.failureReason)
        info.failureReason = '一步更新得到 NaN/Inf';
    end
end

info.uRaw = uRaw;
info.ampInfo = ampInfo;
info.gaugeInfo = gaugeInfo;
info.postGaugeInfo = postGaugeInfo;
info.variant = variant;
end

function v = get_num(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
end
