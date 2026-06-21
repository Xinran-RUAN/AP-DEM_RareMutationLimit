function result = run_wkb_dual_1d(par)
%RUN_WKB_DUAL_1D 实验性双 theta 网格 WKB 求解器。
%
% u(theta,t) 在细网格 theta_u 上推进，W(x,theta,t) 在粗网格 theta_W 上推进。
% 目标：提高相位峰位/曲率分辨率，同时避免把二维变量 W 也加密。
%
% 当前实现范围：
%   1) amplitudeVariant = ''split''：旧版双网格 frozen-H split 更新；
%   2) amplitudeVariant = ''full-eigen-relax''：仍然使用完整的 W 更新格式，
%      即同时包含 x 扩散、theta 扩散和 theta 输运项；唯一变化是 W 粗网格
%      上的 H 直接由同一粗网格的主特征值问题得到，使 eps 很小时 x 方向
%      主导快松弛自然指向对应特征向量。该分支不是把格式近似成
%      逐 theta 的特征向量投影。
% 输出为了兼容后处理，会额外保存 W_on_u 作为变量 W，同时保存 Wcoarse/gridW/opW。

par = local_fill_dual_defaults(par);
if ~strcmpi(local_get_string(par,'timeIntegrator','frozen'), 'frozen')
    error('run_wkb_dual_1d:UnsupportedTimeIntegrator', '双网格版本当前只支持 timeIntegrator=''frozen''。');
end
variant = local_normalize_variant(local_get_string(par,'amplitudeVariant','split'));
if ~any(strcmp(variant, {'split','full-eigen-relax'}))
    error('run_wkb_dual_1d:UnsupportedVariant', ...
        '双网格版本当前只支持 amplitudeVariant=''split'' 或 ''full-eigen-relax''。');
end

% ------------------------- 网格和初值 -------------------------
parW = par;
parW.Ntheta = par.Ntheta_W;
parW.initialGaugeMode = 'none';   % W 粗网格只取振幅初值，统一由细 u 做一次 gauge
[Wc, ~, D_w, Kx, gridW] = model.initial_wkb_1d(parW);
opW = src.build_operators_1d(gridW);

parU = par;
parU.Ntheta = par.Ntheta_u;
gridU = model.make_grid_1d(parU);
opU = src.build_operators_1d(gridU);
D_u = par.D_fun(gridU.theta);

if ~isfield(par, 'u0_fun') || isempty(par.u0_fun)
    u = zeros(1, par.Ntheta_u);
else
    u = par.u0_fun(gridU.theta);
end
u = real(u(:).');
if numel(u) ~= par.Ntheta_u
    error('run_wkb_dual_1d:BadU0', 'u0_fun 返回长度 %d，但 Ntheta_u=%d。', numel(u), par.Ntheta_u);
end
[Wc, u] = src.normalize_gauge(Wc, u, par.eps, 'max');

par.outdir = utils.ensure_dir(par.outdir);
tag = utils.build_tag(par);
prog = utils.progress_init('WKB-DUAL-1D', par);
saveCtl = utils.intermediate_save_init(par, tag, 'wkb');
history = [];
historyEverySteps = local_history_every_steps(par);

% 初始 rho 和质量状态
[rho, nU, rhoInfo, WU] = src.reconstruct_rho_dual_1d(Wc, u, par.eps, opW, opU, par.rhoReconstruction, par.dualWtoUInterp);
massState = src.mass_balance_stats_1d(rho, Kx, opU);
massInfo = local_make_initial_mass_info(massState, par);

% 初始保存
[saveDue0, ~, ~] = utils.intermediate_save_due(saveCtl, 0, 0);
if saveDue0
    H_u_raw = local_solve_H_on_u(D_u, Kx, rho, opU, par, u);
    [H_u, hGaugeInfo] = src.shift_effective_hamiltonian_1d(H_u_raw, par);
    dg = local_compute_dual_diagnostics(Wc, WU, u, D_u, Kx, rho, H_u, nU, gridW, gridU, opU, par, massInfo, massState, rhoInfo, hGaugeInfo);
    vars0 = local_make_save_vars(par, tag, 0, 0, WU, Wc, u, nU, rho, H_u, D_u, D_w, Kx, gridU, gridW, opU, opW, Inf, history, dg, struct(), struct(), variant, massInfo);
    [saveCtl, ~] = utils.intermediate_save_step(saveCtl, par, vars0);
end

% ------------------------- 主循环 -------------------------
t = 0;
step = 0;
residual = Inf;
failed = false;
failureReason = '';
lastInfo = struct();
% 避免浮点舍入导致 T=1000*dt 时多走一个极小的第 1001 步。
timeTol = max(100*eps(max(1,abs(par.T))), 1e-13*max(1,abs(par.T)));

while t < par.T - timeTol && step < par.maxSteps
    remainT = par.T - t;
    if remainT <= timeTol
        t = par.T;
        break;
    end
    dtStep = min(par.dt, remainT);
    if dtStep <= timeTol && remainT <= 10*timeTol
        t = par.T;
        break;
    end
    if dtStep <= 0 || ~isfinite(dtStep)
        failed = true;
        failureReason = 'dtStep 非正或非有限';
        break;
    end

    uOld = u;
    WOld = Wc;

    try
        [rho, ~, rhoInfo0] = src.reconstruct_rho_dual_1d(Wc, u, par.eps, opW, opU, par.rhoReconstruction, par.dualWtoUInterp);
        H_u_raw = local_solve_H_on_u(D_u, Kx, rho, opU, par, u);
        [H_u, hGaugeInfo] = src.shift_effective_hamiltonian_1d(H_u_raw, par);

        uRaw = src.step_phase_1d(uOld, H_u, par.eps, dtStep, opU, par.phaseHamiltonian, par.lfAlpha);
        uRawW = src.interp_periodic_theta_1d(opU.theta, uRaw, opW.theta, par.dualUtoWInterp);

        if strcmp(variant, 'full-eigen-relax')
            % W 更新仍然使用完整 W 方程格式：x 扩散 + theta 扩散 + theta 输运同时隐式/半隐式处理。
            % 为了保持小 eps 时的特征向量松弛结构，W 粗网格上的 H 直接在 W 粗网格求解，
            % 不从 u 细网格插值过来。默认 hamiltonianGauge='none' 时 H_w_update 即 raw H。
            H_w_raw = local_solve_H_on_w(D_w, Kx, rho, opW, par);
            hMode = lower(strrep(local_get_string(par, 'fullEigenRelaxHMode', 'raw'), '_', '-'));
            switch hMode
                case {'raw','eigen-raw','principal-raw'}
                    H_w_update = H_w_raw;
                case {'consistent-gauge','shifted','gauge-consistent'}
                    H_w_update = H_w_raw - local_get_field(hGaugeInfo, 'hamiltonianGaugeShift', 0);
                otherwise
                    error('run_wkb_dual_1d:BadFullEigenRelaxHMode', ...
                        '未知 fullEigenRelaxHMode=%s。', hMode);
            end
            [WStep, ampInfo] = src.update_w_full_eigen_relax_fd_1d( ...
                WOld, uRawW, D_w, Kx, rho, H_w_update, H_w_raw, par.eps, dtStep, opW, par);
        else
            H_w_raw = src.interp_periodic_theta_1d(opU.theta, H_u_raw, opW.theta, par.dualHtoWInterp);
            H_w = H_w_raw - local_get_field(hGaugeInfo, 'hamiltonianGaugeShift', 0);
            [WStep, ampInfo] = src.update_w_split_fd_1d(WOld, uRawW, D_w, Kx, rho, H_w, par.eps, dtStep, opW, par);
        end

        [WTrial, uTrial] = src.normalize_gauge(WStep, uRaw, par.eps, 'max');

        massInfoTrial = massInfo;
        massStateTrial = massState;
        if local_mass_correction_requested(par)
            [rhoTrialForMass, ~] = src.reconstruct_rho_dual_1d(WTrial, uTrial, par.eps, opW, opU, par.rhoReconstruction, par.dualWtoUInterp);
            [massScale, massInfoTrial] = src.mass_balance_correction_1d( ...
                rhoTrialForMass, Kx, opU, par.eps, dtStep, ...
                massState.mass, massState.reactionIntegral, par);
            if massInfoTrial.applied
                WTrial = massScale * WTrial;
                rhoTrialForMass = massScale * rhoTrialForMass;
            end
            massStateTrial = src.mass_balance_stats_1d(rhoTrialForMass, Kx, opU);
            massInfoTrial = local_finalize_mass_info(massInfoTrial, massStateTrial);
        end

        if any(~isfinite(WTrial(:))) || any(~isfinite(uTrial(:)))
            error('trial W/u 出现 NaN 或 Inf。');
        end

        % 接受本步
        Wc = WTrial;
        u = uTrial;
        massInfo = massInfoTrial;
        massState = massStateTrial;
        t = t + dtStep;
        if abs(t - par.T) <= timeTol
            t = par.T;
        end
        step = step + 1;
        residual = max(abs(u(:) - uOld(:))) / max(dtStep, realmin);

        % residual<=tol 时，若 stopByResidual=false，不再每一步都强制记录/打印；
        % 否则收敛后会出现大量逐步输出，看起来像后期异常。
        hitTolForHistory = local_get_bool(par, 'stopByResidual', false) && residual <= par.tol;
        historyDue = (step == 1 || mod(step, historyEverySteps) == 0 || hitTolForHistory);
        [saveDue, ~, ~] = utils.intermediate_save_due(saveCtl, t, step);

        if historyDue || saveDue
            [rhoNow, nNow, rhoInfoNow, WU] = src.reconstruct_rho_dual_1d(Wc, u, par.eps, opW, opU, par.rhoReconstruction, par.dualWtoUInterp);
            HNowRaw = local_solve_H_on_u(D_u, Kx, rhoNow, opU, par, u);
            [HNow, hGaugeInfoNow] = src.shift_effective_hamiltonian_1d(HNowRaw, par);
            massState = src.mass_balance_stats_1d(rhoNow, Kx, opU);
            massInfo = local_finalize_mass_info(massInfo, massState);
            dg = local_compute_dual_diagnostics(Wc, WU, u, D_u, Kx, rhoNow, HNow, nNow, gridW, gridU, opU, par, massInfo, massState, rhoInfoNow, hGaugeInfoNow);
            dg = local_merge_structs(dg, ampInfo);
            if historyDue
                history = utils.append_history(history, t, residual, dg, struct(), struct('dt',dtStep));
            end
            lastInfo = local_progress_info(dg);
            if saveDue
                vars = local_make_save_vars(par, tag, t, step, WU, Wc, u, nNow, rhoNow, HNow, D_u, D_w, Kx, gridU, gridW, opU, opW, residual, history, dg, struct(), struct('dt',dtStep), variant, massInfo);
                [saveCtl, ~] = utils.intermediate_save_step(saveCtl, par, vars);
            end
        end

        prog = utils.progress_update(prog, step, t, residual, lastInfo);
        if par.stopByResidual && isfinite(residual) && residual <= par.tol
            break;
        end
    catch ME
        failed = true;
        failureReason = local_error_report(ME);
        if par.verbose
            fprintf(2, '\n[WKB-DUAL-1D] 失败：%s\n', failureReason);
        end
        break;
    end
end

% ------------------------- 最终输出 -------------------------
try
    [rho, nU, rhoInfo, WU] = src.reconstruct_rho_dual_1d(Wc, u, par.eps, opW, opU, par.rhoReconstruction, par.dualWtoUInterp);
    Hraw = local_solve_H_on_u(D_u, Kx, rho, opU, par, u);
    [H, hGaugeInfoFinal] = src.shift_effective_hamiltonian_1d(Hraw, par);
    massState = src.mass_balance_stats_1d(rho, Kx, opU);
    massInfo = local_finalize_mass_info(massInfo, massState);
    diagInfo = local_compute_dual_diagnostics(Wc, WU, u, D_u, Kx, rho, H, nU, gridW, gridU, opU, par, massInfo, massState, rhoInfo, hGaugeInfoFinal);
catch ME
    WU = NaN(opU.nx, opU.Ntheta); nU = WU; rho = NaN(opU.nx,1); H = NaN(1,opU.Ntheta);
    diagInfo = struct('failureReason', ['final diagnostics failed: ' ME.message]);
    residual = Inf;
end

diagInfo.failed = failed;
diagInfo.failureReason = failureReason;

result = struct('W', WU, 'Wcoarse', Wc, 'u', u, 'n', nU, 'rho', rho, 'H', H, ...
                'D', D_u, 'Dcoarse', D_w, 'Kx', Kx, 'grid', gridU, 'gridW', gridW, ...
                'op', opU, 'opW', opW, 'par', par, 't', t, 'step', step, 'residual', residual, ...
                'history', history, 'diagnostics', diagInfo, 'tag', tag, 'variant', variant, ...
                'massInfo', massInfo, 'dualThetaInfo', local_dual_info(par, gridU, gridW));

try
    finalVars = local_make_save_vars(par, tag, t, step, WU, Wc, u, nU, rho, H, D_u, D_w, Kx, gridU, gridW, opU, opW, residual, history, diagInfo, struct(), struct(), variant, massInfo);
    [saveCtl, ~] = utils.intermediate_save_force_latest(saveCtl, par, finalVars, 'final-or-failed-checkpoint'); %#ok<NASGU>
catch ME
    warning('run_wkb_dual_1d:FinalCheckpointFailed', '最终 checkpoint 保存失败：%s', ME.message);
end

utils.progress_finish(prog, step, t, residual, local_progress_info(diagInfo));
if local_get_bool(par, 'saveFinalResult', true)
    outfile = utils.solver_file(par.outdir, 'wkb', 'final');
    save(outfile, '-struct', 'result');
    if par.verbose
        fprintf('[WKB-DUAL-1D] 结果已保存: %s\n', outfile);
    end
end
end

% =====================================================================
% 局部函数
% =====================================================================

function par = local_fill_dual_defaults(par)
% 双网格分支不再在这里设置隐藏默认值。所有相关参数必须由
% run_wkb_main_1d.m 显式给出，经 run_case_1d 传入。
required = {'Ntheta_W','Ntheta_u','dualWtoUInterp','dualUtoWInterp','dualHtoWInterp', ...
    'HSolverMode','HInterpND','HLocalRadius','dualHSolverModeW','fullEigenRelaxHMode', ...
    'eigenSolver','maxSteps','tol','verbose','stopByResidual','lfAlpha'};
missing = {};
for i = 1:numel(required)
    f = required{i};
    if ~isfield(par, f) || isempty(par.(f))
        if strcmp(f, 'Ntheta_W') && isfield(par, 'Ntheta') && ~isempty(par.Ntheta)
            % run_case_1d 已约定 par.Ntheta 是 W 粗网格点数；Ntheta_W=[] 时在这里明示等于 par.Ntheta。
            continue;
        end
        missing{end+1} = f; %#ok<AGROW>
    end
end
if ~isempty(missing)
    error('run_wkb_dual_1d:MissingExplicitDualParameters', ...
        '双 theta 网格缺少显式参数：%s。请在 run_wkb_main_1d.m 中写明。', strjoin(missing, ', '));
end

if ~isfield(par, 'Ntheta_W') || isempty(par.Ntheta_W)
    par.Ntheta_W = par.Ntheta;
end
par.Ntheta_W = round(par.Ntheta_W);
par.Ntheta_u = round(par.Ntheta_u);
par.Ntheta = par.Ntheta_W;
par.useDualThetaGrid = true;
end

function H = local_solve_H_on_u(Du, Kx, rho, opU, par, u)
parH = par;
mode = lower(strrep(local_get_string(par, 'HSolverMode', 'direct-unique-D'), '_', '-'));
if any(strcmp(mode, {'d-pchip-local-correct','dinterp-local-correct'}))
    [~, kPeak] = max(real(u(:)));
    r = max(0, round(local_get_field(par, 'HLocalRadius', 4)));
    K = numel(u);
    idx = mod((kPeak-r:kPeak+r)-1, K) + 1;
    parH.HLocalCorrectionIndices = unique(idx);
end
H = src.solve_H_fd_1d(Du, Kx, rho, opU, parH);
end

function H = local_solve_H_on_w(Dw, Kx, rho, opW, par)
% W 粗网格的 H。full-eigen-relax 分支要求 H 与 W 粗网格上的 x-算子一致，
% 因此默认直接在粗网格 theta_W 上求，而不是从 u 细网格插值。
parH = par;
parH.HSolverMode = local_get_string(par, 'dualHSolverModeW', 'direct-unique-D');
H = src.solve_H_fd_1d(Dw, Kx, rho, opW, parH);
end

function dg = local_compute_dual_diagnostics(Wc, WU, u, D, Kx, rho, H, nU, gridW, gridU, opU, par, massInfo, massState, rhoInfo, hInfo)
dg = struct();
dg.dualThetaGrid = true;
dg.Ntheta_u = gridU.Ntheta;
dg.Ntheta_W = gridW.Ntheta;
dg.theta_u = gridU.theta;
dg.theta_W = gridW.theta;
dg.rho_max = max(rho(:));
dg.rho_min = min(rho(:));
dg.Hmin = min(H(:));
dg.Hmax = max(H(:));
P = trapz(opU.x(:), nU, 1).';
dg.P = real(P(:));
[~, ip] = max(dg.P);
[~, iu] = max(u(:));
dg.theta_direct = gridU.theta(ip);
dg.theta_wkb = gridU.theta(iu);
if isfield(par, 'theta_m')
    dg.err_wkb_to_m = utils.periodic_distance(dg.theta_wkb, par.theta_m);
    dg.err_direct_to_m = utils.periodic_distance(dg.theta_direct, par.theta_m);
end
massP = opU.dtheta * sum(max(dg.P,0));
if massP > 0
    dist = arrayfun(@(th) utils.periodic_distance(th, dg.theta_wkb), gridU.theta(:));
    dg.sigma_theta = sqrt(sum((dist.^2).*max(dg.P(:),0))*opU.dtheta / massP);
else
    dg.sigma_theta = NaN;
end
dg.WcoarseMax = max(Wc(:));
dg.WcoarseMin = min(Wc(:));
dg.WonUMax = max(WU(:));
dg.WonUMin = min(WU(:));
if isstruct(rhoInfo)
    dg.rhoReconstructLogMax = local_get_field(rhoInfo, 'logrhoMax', NaN);
end
if isstruct(hInfo)
    dg = local_merge_structs(dg, hInfo);
end
dg = local_attach_mass_info(dg, massInfo, massState);
end

function vars = local_make_save_vars(par, tag, t, step, WU, Wc, u, n, rho, H, D, Dw, Kx, grid, gridW, op, opW, residual, history, dg, resInfo, dtInfo, variant, massInfo)
vars = struct();
vars.W = WU;                 % 兼容旧后处理：W 已插值到 u 细网格
vars.Wcoarse = Wc;
vars.u = u;
vars.n = n;
vars.rho = rho;
vars.H = H;
vars.D = D;
vars.Dcoarse = Dw;
vars.Kx = Kx;
vars.grid = grid;
vars.gridW = gridW;
vars.op = op;
vars.opW = opW;
vars.par = par;
vars.tag = tag;
vars.t = t;
vars.step = step;
vars.residual = residual;
vars.history = history;
vars.diagnostics = dg;
vars.resInfo = resInfo;
vars.dtInfo = dtInfo;
vars.variant = variant;
vars.dualThetaInfo = local_dual_info(par, grid, gridW);
if nargin >= 24 && isstruct(massInfo)
    vars.massInfo = massInfo;
end
vars.createdBy = 'run_wkb_dual_1d intermediate output';
end

function info = local_dual_info(par, gridU, gridW)
info = struct();
info.enabled = true;
info.Ntheta_u = gridU.Ntheta;
info.Ntheta_W = gridW.Ntheta;
info.theta_u = gridU.theta;
info.theta_W = gridW.theta;
info.WtoUInterp = local_get_string(par, 'dualWtoUInterp', 'pchip');
info.UtoWInterp = local_get_string(par, 'dualUtoWInterp', 'pchip');
info.HtoWInterp = local_get_string(par, 'dualHtoWInterp', 'pchip');
info.HSolverMode = local_get_string(par, 'HSolverMode', '');
info.dualHSolverModeW = local_get_string(par, 'dualHSolverModeW', '');
info.fullEigenRelaxHMode = local_get_string(par, 'fullEigenRelaxHMode', '');
end

function tf = local_mass_correction_requested(par)
tf = local_get_bool(par, 'massCorrectionDuringSolve', false);
if isstruct(par) && isfield(par, 'massCorrection') && ~isempty(par.massCorrection)
    mode = lower(char(par.massCorrection));
    if any(strcmp(mode, {'none','off','false','0','no'})), tf = false; end
    if any(strcmp(mode, {'on','true','1','yes','balance','balance-be','balance-trapezoid'})), tf = true; end
end
end

function info = local_make_initial_mass_info(massState, par)
info = struct();
info.enabled = local_mass_correction_requested(par);
info.formula = local_get_string(par, 'massCorrectionFormula', 'trapezoid');
info.applied = false;
info.scale = 1;
info.reason = 'initial';
info.correctedMass = massState.mass;
info.correctedReactionIntegral = massState.reactionIntegral;
info.massAfter = massState.mass;
info.reactionIntegralAfter = massState.reactionIntegral;
info.trialMass = massState.mass;
end

function info = local_finalize_mass_info(info, massState)
if ~isstruct(info), info = struct(); end
info.massAfter = massState.mass;
info.reactionIntegralAfter = massState.reactionIntegral;
info.correctedMass = massState.mass;
info.correctedReactionIntegral = massState.reactionIntegral;
info.rhoMaxAfter = massState.rhoMax;
info.rhoMinAfter = massState.rhoMin;
end

function diag = local_attach_mass_info(diag, massInfo, massState)
diag.mass_total = massState.mass;
diag.mass_balance_rhs = massState.reactionIntegral;
if isstruct(massInfo)
    diag.mass_correction_enabled = local_get_bool(massInfo, 'enabled', false);
    diag.mass_correction_applied = local_get_bool(massInfo, 'applied', false);
    diag.mass_correction_scale = local_get_field(massInfo, 'scale', NaN);
    diag.mass_trial = local_get_field(massInfo, 'trialMass', NaN);
    diag.mass_corrected = local_get_field(massInfo, 'correctedMass', massState.mass);
end
end

function historyEverySteps = local_history_every_steps(par)
if isfield(par, 'historyEverySteps') && ~isempty(par.historyEverySteps) && isnumeric(par.historyEverySteps) && isscalar(par.historyEverySteps) && par.historyEverySteps > 0
    historyEverySteps = max(1, round(par.historyEverySteps));
elseif isfield(par, 'historyEveryTime') && ~isempty(par.historyEveryTime) && isnumeric(par.historyEveryTime) && isscalar(par.historyEveryTime) && par.historyEveryTime > 0
    historyEverySteps = max(1, round(par.historyEveryTime / par.dt));
else
    historyEverySteps = max(1, round(0.1 / par.dt));
end
end

function variant = local_normalize_variant(name)
name = lower(char(name));
switch name
    case {'split','split-wkb'}
        variant = 'split';
    case {'full-eigen-relax','full_eigen_relax','fulleigenrelax','eigen-relax-full','eigen_relax_full','dual-full-eigen-relax'}
        variant = 'full-eigen-relax';
    otherwise
        variant = name;
end
end

function info = local_progress_info(dg)
info = struct();
if isempty(dg) || ~isstruct(dg), return; end
info.theta_wkb = local_get_field(dg,'theta_wkb',NaN);
info.theta_direct = local_get_field(dg,'theta_direct',NaN);
info.err_wkb_to_m = local_get_field(dg,'err_wkb_to_m',NaN);
info.sigma_theta = local_get_field(dg,'sigma_theta',NaN);
info.rho_max = local_get_field(dg,'rho_max',NaN);
info.mass_total = local_get_field(dg,'mass_total',NaN);
info.mass_correction_scale = local_get_field(dg,'mass_correction_scale',NaN);
info.Hmin = local_get_field(dg,'Hmin',NaN);
info.xEigenResidualInfMax = local_get_field(dg,'xEigenResidualInfMax',NaN);
info.solveResidualInf = local_get_field(dg,'solveResidualInf',NaN);
end


function msg = local_error_report(ME)
msg = ME.message;
try
    if ~isempty(ME.identifier)
        msg = sprintf('%s [%s]', msg, ME.identifier);
    end
    if ~isempty(ME.stack)
        st = ME.stack(1);
        msg = sprintf('%s at %s:%d', msg, st.name, st.line);
    end
catch
    msg = ME.message;
end
end

function out = local_merge_structs(varargin)
out = struct();
for i = 1:nargin
    s = varargin{i};
    if ~isstruct(s), continue; end
    f = fieldnames(s);
    for j = 1:numel(f)
        out.(f{j}) = s.(f{j});
    end
end
end

function v = local_get_field(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
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
