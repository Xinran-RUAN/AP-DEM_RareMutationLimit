function result = run_wkb_timeap_1d(par)
%RUN_WKB_TIMEAP_1D  v11 实验版 WKB 时间积分器。
%
%  本函数实现用户要求测试的两个方向：
%
%  方案 1：par.timeIntegrator = 'implicit-coupled'
%    全隐式耦合 Picard 时间步。每个大时间步内用 Picard 迭代同时更新
%    rho^{n+1}, H^{n+1}, u^{n+1}, W^{n+1}。这样做的目的不是求稳态，
%    而是在 time evolution 中避免简单冻结 H^n 带来的 O(dt) 相位误差。
%
%  方案 3：par.timeIntegrator = 'projective'
%    projective / micro-macro 时间积分。先用若干 micro steps 解析 O(eps)
%    快松弛层，然后对慢变量做宏步外推。默认只外推 phase u 的形状，
%    W 作为快变量保留在 micro-relaxed 状态，避免 log(W) 外推导致爆炸。
%
%  两个方案都保持 WKB 框架：主变量始终是 W 和 u。代码中出现的重构密度
%  rho 只用于计算非局部项和诊断，不作为主要求解变量。
%
%  重要提醒：这两个算法是研究/排错版，不是已经证明的 time-AP 格式。
%  它们的作用是帮助判断：小 epsilon 下的失败到底来自 H(t) 冻结误差、
%  快层没有解析，还是来自 rho 重构/相位峰欠解析。

par = local_preconfigure(par);
[W, u, D, Kx, grid] = model.initial_wkb_1d(par);
op = src.build_operators_1d(grid);
variant = local_normalize_variant(par.amplitudeVariant);
timeIntegrator = lower(char(par.timeIntegrator));

par.outdir = utils.ensure_dir(par.outdir);
tag = utils.build_tag(par);
prog = utils.progress_init('WKB-1D', par);
lp = utils.liveplot_init_1d(par, 'wkb');
livePlotEverySteps = local_every_steps(par, 'livePlotEverySteps', 'livePlotEveryTime', 1.0);
historyEverySteps = local_every_steps(par, 'historyEverySteps', 'historyEveryTime', 0.5);

history = [];
% 中间输出控制器。实验版 time-AP/projective 计算通常很慢，
% 如果只在最终时刻保存，一旦中途停止就无法检查前面的计算。
% 因此这里维护 full snapshot、latest checkpoint 和 history checkpoint 三类输出。
saveCtl = local_intermediate_save_init(par, tag);

t = 0;
step = 0;
residual = Inf;
failed = false;
failureReason = '';
lastInfo = struct();
lastPlotState = struct();

[rhoInitialForMass, ~] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
massState = src.mass_balance_stats_1d(rhoInitialForMass, Kx, op);
massInfoCurrent = local_make_initial_mass_info(massState, par);

% 可选保存初始时刻 t=0。time-AP 路径原来只在第一个宏步后保存，
% 若 snapshotTimes 含 0 或 saveInitialSnapshot=true，后处理会缺少初值。
% 这里仅写文件，不改变 W/u/rho/H 或时间积分格式。
try
    [rho0, n0] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
    H0raw = src.solve_H_fd_1d(D, Kx, rho0, op);
    [H0, hGaugeInfo0] = src.shift_effective_hamiltonian_1d(H0raw, par);
    dg0 = src.compute_diagnostics_1d(W, u, D, Kx, rho0, H0, par.eps, op, par, variant);
    dg0 = local_merge_structs(dg0, hGaugeInfo0);
    dg0 = local_attach_mass_info(dg0, massInfoCurrent, massState);
    [saveCtl, didSaveInitial] = local_maybe_save_intermediate(saveCtl, par, tag, ...
        t, step, W, u, n0, rho0, H0, D, Kx, grid, op, residual, history, dg0, struct(), struct(), struct('saveReason','initial'), massInfoCurrent);
    if didSaveInitial && par.verbose && local_get_bool(par,'printIntermediateSaveMessage',false)
        fprintf('[WKB-1D] 已保存初始中间结果: step=0, t=0\n');
    end
catch ME
    warning('run_wkb_timeap_1d:InitialSnapshotFailed', '初始 time-AP/WKB 快照保存失败：%s', ME.message);
end

while t < par.T && step < par.maxSteps
    Wold = W;
    uold = u;
    [rhoOld, ~] = src.reconstruct_rho_1d(Wold, uold, par.eps, op, par.rhoReconstruction);
    HrawOld = src.solve_H_fd_1d(D, Kx, rhoOld, op);
    [Hold, ~] = src.shift_effective_hamiltonian_1d(HrawOld, par);
    [dtStep, dtInfo] = src.select_wkb_time_step_1d(par, t, Wold, uold, D, Kx, rhoOld, Hold, op, variant);
    if dtStep <= 0 || ~isfinite(dtStep)
        failed = true;
        failureReason = '时间步 dtStep 非正或非有限';
        break;
    end

    switch timeIntegrator
        case {'implicit-coupled','fully-implicit-coupled','coupled-implicit'}
            [Wtry, utry, stepInfo] = local_step_implicit_coupled(Wold, uold, D, Kx, dtStep, op, par, variant);
        case {'projective','micro-macro','projective-micro'}
            [Wtry, utry, stepInfo] = local_step_projective(Wold, uold, D, Kx, dtStep, op, par, variant);
        otherwise
            error('run_wkb_timeap_1d should not be called for timeIntegrator="%s".', timeIntegrator);
    end

    if stepInfo.failed
        failed = true;
        failureReason = stepInfo.failureReason;
        break;
    end

    % 每个宏步结束后可选进行质量平衡投影。这里只缩放 W，保持 u 不变。
    massInfoStep = massInfoCurrent;
    massStateStep = massState;
    if local_mass_correction_requested(par)
        try
            [rhoTryMass, ~] = src.reconstruct_rho_1d(Wtry, utry, par.eps, op, par.rhoReconstruction);
            [massScale, massInfoStep] = src.mass_balance_correction_1d( ...
                rhoTryMass, Kx, op, par.eps, dtStep, ...
                massState.mass, massState.reactionIntegral, par);
            if massInfoStep.applied
                Wtry = massScale * Wtry;
                rhoTryMass = massScale * rhoTryMass;
            end
            massStateStep = src.mass_balance_stats_1d(rhoTryMass, Kx, op);
            massInfoStep = local_finalize_mass_info(massInfoStep, massStateStep);
        catch ME
            failed = true;
            failureReason = ['质量平衡修正失败: ' ME.message];
            break;
        end
    end

    if ~local_mass_correction_requested(par)
        [rhoTryMass, ~] = src.reconstruct_rho_1d(Wtry, utry, par.eps, op, par.rhoReconstruction);
        massStateStep = src.mass_balance_stats_1d(rhoTryMass, Kx, op);
        massInfoStep = local_finalize_mass_info(massInfoStep, massStateStep);
    end

    W = Wtry;
    u = utry;
    massState = massStateStep;
    massInfoCurrent = massInfoStep;
    t = t + dtStep;
    step = step + 1;

    % 每个宏步结束后，用当前状态重新计算 rho/H/diagnostics。
    [rhoNow, nNow, rhoInfo] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
    HrawNow = src.solve_H_fd_1d(D, Kx, rhoNow, op);
    [HNow, hGaugeInfo] = src.shift_effective_hamiltonian_1d(HrawNow, par);

    residualMode = lower(local_get_string(par, 'residualMode', 'wkb-steady'));
    [updateResidual, resInfo] = src.compute_wkb_residual_1d(W, u, Wold, uold, par.eps, dtStep, op, 'scaled-gauge-invariant');
    steadyInfo = src.compute_wkb_steady_residual_1d(W, u, D, Kx, rhoNow, HNow, par.eps, op, par, variant);
    if any(strcmp(residualMode, {'wkb-steady','steady','steady-wkb'}))
        residual = steadyInfo.res_wkb_steady;
    else
        residual = updateResidual;
    end

    dtInfo.dt = dtStep;
    dtInfo.timeIntegratorCode = local_integrator_code(timeIntegrator);
    if isfield(stepInfo,'numRetries'), dtInfo.numRetries = stepInfo.numRetries; end

    dg = src.compute_diagnostics_1d(W, u, D, Kx, rhoNow, HNow, par.eps, op, par, variant);
    monInfo = src.monitor_wkb_small_eps_1d(W, u, Wold, uold, par.eps, dtStep, op, par, stepInfo.ampInfo, dtInfo, steadyInfo);
    % 字段名兼容 utils.append_history / progress_update。
    monInfo.neighborJumpOverEps = monInfo.maxNeighborJumpOverEps;
    if isfield(monInfo,'WMax'), monInfo.W_max = monInfo.WMax; end
    if isfield(stepInfo,'ampInfo') && isstruct(stepInfo.ampInfo)
        ampInfo = stepInfo.ampInfo;
    else
        ampInfo = struct();
    end
    massState = src.mass_balance_stats_1d(rhoNow, Kx, op);
    massInfoCurrent = local_finalize_mass_info(massInfoCurrent, massState);
    dg = local_merge_structs(dg, hGaugeInfo, steadyInfo, monInfo, ampInfo, stepInfo);
    dg = local_attach_mass_info(dg, massInfoCurrent, massState);
    dg.rhoReconstructLogMax = get_field(rhoInfo,'logrhoMax',NaN);
    dg.rhoReconstructNumCappedHigh = get_field(rhoInfo,'numCappedHigh',get_field(rhoInfo,'numCappedRhoHigh',0));

    if step == 1 || mod(step, historyEverySteps) == 0 || residual <= par.tol || get_field(stepInfo,'forceHistory',0) > 0
        history = utils.append_history(history, t, residual, dg, resInfo, dtInfo);
        lastInfo = local_progress_info(dg, resInfo, dtInfo, timeIntegrator);
        lastPlotState = local_make_plot_state(W, u, nNow, rhoNow, HNow, Kx, op, par, t, step, residual, dg);
    end

    % 已接受的宏步结束后保存中间结果。这里保存的是当前 WKB 状态，
    % 不保存被拒绝的 trial step，也不改变数值格式。
    [saveCtl, didSaveIntermediate] = local_maybe_save_intermediate(saveCtl, par, tag, ...
        t, step, W, u, nNow, rhoNow, HNow, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, stepInfo, massInfoCurrent);
    if didSaveIntermediate && par.verbose && local_get_bool(par,'printIntermediateSaveMessage',false)
        fprintf('[WKB-1D] 已保存中间结果: step=%d, t=%.6g\n', step, t);
    end

    prog = utils.progress_update(prog, step, t, residual, lastInfo);

    if lp.enabled && ~isempty(fieldnames(lastPlotState)) && ...
            (step == 1 || mod(step, livePlotEverySteps) == 0 || residual <= par.tol)
        lastPlotState.saveFrame = isfield(par, 'livePlotSave') && par.livePlotSave;
        lp = utils.liveplot_update_1d(lp, 'wkb', lastPlotState);
        if isfield(par, 'livePlotPause') && par.livePlotPause > 0
            pause(par.livePlotPause);
        end
    end

    if par.stopByResidual && isfinite(residual) && residual <= par.tol
        break;
    end

    if any(~isfinite(W(:))) || any(~isfinite(u(:))) || ~isfinite(residual)
        failed = true;
        failureReason = '宏步结束后 W/u/residual 出现 NaN/Inf';
        break;
    end
end

try
    [rho, n] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
    Hraw = src.solve_H_fd_1d(D, Kx, rho, op);
    [H, hGaugeInfoFinal] = src.shift_effective_hamiltonian_1d(Hraw, par);
    massState = src.mass_balance_stats_1d(rho, Kx, op);
    massInfoCurrent = local_finalize_mass_info(massInfoCurrent, massState);
    diagInfo = src.compute_diagnostics_1d(W, u, D, Kx, rho, H, par.eps, op, par, variant);
    steadyFinal = src.compute_wkb_steady_residual_1d(W, u, D, Kx, rho, H, par.eps, op, par, variant);
    diagInfo = local_merge_structs(diagInfo, hGaugeInfoFinal, steadyFinal);
    diagInfo = local_attach_mass_info(diagInfo, massInfoCurrent, massState);
catch ME
    rho = NaN(op.nx,1); n = NaN(op.nx,op.Ntheta); H = NaN(1,op.Ntheta);
    diagInfo = struct('theta_wkb',NaN,'err_wkb_to_m',NaN,'sigma_theta',NaN, ...
        'rho_max',Inf,'gammaR',Inf,'failureReason',['final diagnostics failed: ' ME.message]);
    residual = Inf;
end
if failed
    diagInfo.failed = true;
    diagInfo.failureReason = failureReason;
    if par.verbose
        fprintf('\n[WKB-1D] 警告：实验版时间积分器提前停止。原因: %s\n', failureReason);
    end
else
    diagInfo.failed = false;
    diagInfo.failureReason = '';
end


% 无论正常结束、达到终止条件、失败退出，最后都补写一次 latest checkpoint。
% 这样即使没有跑到 par.T，也能在输出目录中找到最近状态。
try
    [saveCtl, ~] = local_force_save_latest(saveCtl, par, tag, ...
        t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, diagInfo, struct(), struct(), struct(), massInfoCurrent); %#ok<NASGU>
catch ME
    warning('run_wkb_timeap_1d:FinalCheckpointFailed', ...
        '最终 checkpoint 保存失败：%s', ME.message);
end

result = struct('W',W,'u',u,'n',n,'rho',rho,'H',H,'D',D,'Kx',Kx,'grid',grid, ...
                'op',op,'par',par,'t',t,'step',step,'residual',residual, ...
                'history',history,'diagnostics',diagInfo,'tag',tag,'variant',variant, ...
                'timeIntegrator',timeIntegrator,'massInfo',massInfoCurrent);

utils.progress_finish(prog, step, t, residual, local_progress_info(diagInfo, struct(), struct(), timeIntegrator));
if local_get_bool(par, 'saveFinalResult', true)
    outfile = utils.solver_file(par.outdir, 'wkb', 'final');
    save(outfile, '-struct', 'result');
    if par.verbose
        fprintf('[WKB-1D] 结果已保存: %s\n', outfile);
    end
else
    if par.verbose
        fprintf('[WKB-1D] saveFinalResult=false，仅返回结果并保留已写出的 snapshot。\n');
    end
end
utils.liveplot_finish_1d(lp);
end

%==========================================================================
% 方案 1：全隐式耦合 Picard
%==========================================================================
function [Wnew, unew, info] = local_step_implicit_coupled(Wold, uold, D, Kx, dt, op, par, variant)
% 在一个宏步内求解近似固定点：
%   rho^* = R[W^*,u^*], H^*=H(rho^*),
%   (W^*,u^*) = Step(W^n,u^n; rho^*,H^*).
% 这不是直接求稳态，而是全隐式时间步。Picard 迭代收敛后，H 使用的是
% 当前时间层的 rho，而不是旧时间层的 rho。
info = struct('failed',false,'failureReason','','implicitIter',0, ...
              'implicitRelChange',Inf,'implicitConverged',0,'numRetries',0);
maxIter = round(local_get_num(par,'implicitMaxIter',12));
tol = local_get_num(par,'implicitTol',1e-7);
relax = local_get_num(par,'implicitRelax',0.7);
minIter = round(local_get_num(par,'implicitMinIter',2));
relax = min(max(relax, 0.05), 1.0);

Wguess = Wold;
uguess = uold;
Wnew = Wold;
unew = uold;
lastStepInfo = struct();

for it = 1:maxIter
    [rhoGuess, ~] = src.reconstruct_rho_1d(Wguess, uguess, par.eps, op, par.rhoReconstruction);
    if any(~isfinite(rhoGuess(:)))
        info.failed = true;
        info.failureReason = 'implicit Picard 中 rhoGuess 出现 NaN/Inf';
        return;
    end
    Hraw = src.solve_H_fd_1d(D, Kx, rhoGuess, op);
    if any(~isfinite(Hraw(:)))
        info.failed = true;
        info.failureReason = 'implicit Picard 中 H 特征值计算失败';
        return;
    end
    [Hguess, hInfo] = src.shift_effective_hamiltonian_1d(Hraw, par);
    [Wcand, ucand, oneInfo] = src.wkb_one_step_frozenH_1d(Wold, uold, D, Kx, rhoGuess, Hguess, par.eps, dt, op, par, variant);
    lastStepInfo = oneInfo;
    if oneInfo.failed
        info.failed = true;
        info.failureReason = ['implicit Picard 一步算子失败: ' oneInfo.failureReason];
        info.ampInfo = oneInfo.ampInfo;
        return;
    end

    % Picard 欠松弛。这里的混合仅用于非线性迭代，不是额外物理模型。
    Wmix = (1-relax)*Wguess + relax*Wcand;
    umix = (1-relax)*uguess + relax*ucand;
    Wmix = max(real(Wmix), 0);
    [Wmix, umix] = src.normalize_gauge(Wmix, umix, par.eps, 'max');

    relU = max(abs(umix(:)-uguess(:))) / max(1, max(abs(umix(:))));
    relW = max(abs(Wmix(:)-Wguess(:))) / max(1, max(abs(Wmix(:))));
    rel = max(relU, relW);

    Wguess = Wmix;
    uguess = umix;
    Wnew = Wguess;
    unew = uguess;

    info.implicitIter = it;
    info.implicitRelChange = rel;
    info.HrawMinImplicit = get_field(hInfo,'HrawMin',NaN);
    info.HrawMaxImplicit = get_field(hInfo,'HrawMax',NaN);

    if it >= minIter && rel < tol
        info.implicitConverged = 1;
        break;
    end
end

if info.implicitConverged == 0
    failTol = local_get_num(par,'implicitFailureTol',5e-2);
    if local_get_bool(par,'implicitRetryOnFail',true) && info.implicitRelChange > failTol
        info.failed = true;
        info.failureReason = sprintf('implicit Picard 未充分收敛: rel=%.3g', info.implicitRelChange);
        return;
    end
end

if isfield(lastStepInfo,'ampInfo')
    info.ampInfo = lastStepInfo.ampInfo;
else
    info.ampInfo = struct();
end
info.forceHistory = double(info.implicitConverged == 0);
end

%==========================================================================
% 方案 3：projective / micro-macro
%==========================================================================
function [Wnew, unew, info] = local_step_projective(Wold, uold, D, Kx, dtMacro, op, par, variant)
% 先用 micro steps 解析快层，再对慢变量投影。默认 mode='u-only'：
%   - W 在 micro-relaxed 状态停止，不做大步外推；
%   - u 的形状用最后两个 micro steps 的差分外推；
% 这样比同时外推 log(W) 更稳健，也更符合“W 是快变量、u 是慢相位”的
% micro-macro 思路。
%
% 重要的稳健性处理：本函数内部有多处可能提前 return 的保护分支，
% 例如 micro step 或 corrector step 失败时会设置 info.failed=true。
% MATLAB 要求所有输出变量在 return 前必须已经赋值；否则会出现
%   Output argument "Wnew" not assigned
% 的报错。这里在函数开头先把输出初始化为旧状态。若后续失败，
% 主循环会收到 info.failed 并优雅停止，同时 checkpoint_wkb_latest.mat 中
% 仍然保留失败前的最后一个有效状态。
Wnew = Wold;
unew = uold;

info = struct('failed',false,'failureReason','','projectiveMicroStepsUsed',0, ...
              'projectiveMicroDt',NaN,'projectiveRemainingDt',NaN, ...
              'projectiveModeCode',0,'numRetries',0);

mode = lower(local_get_string(par,'projectiveMode','u-only'));
M = max(2, round(local_get_num(par,'projectiveMicroSteps',8)));
corrector = max(0, round(local_get_num(par,'projectiveCorrectorSteps',0)));

% 根据当前 H/rho 估计 micro dt。这里允许 macro dt 与 eps 无关，但 micro dt
% 必须能解析 O(eps) 快层。注意总的 micro+corrector 时间不能超过宏步 dtMacro，
% 否则 projective 步会悄悄多走物理时间；因此后面会再次截断 dtm。
[rho0, ~] = src.reconstruct_rho_1d(Wold, uold, par.eps, op, par.rhoReconstruction);
Hraw0 = src.solve_H_fd_1d(D, Kx, rho0, op);
[H0, ~] = src.shift_effective_hamiltonian_1d(Hraw0, par);
scale = max([1, max(abs(H0(:))), max(max(Kx(:)-rho0(:),0))]);
dtm = local_get_num(par,'projectiveMicroDtFactor',0.25) * par.eps / scale;
dtm = max(dtm, local_get_num(par,'projectiveMinMicroDt',1e-12));
dtm = min(dtm, local_get_num(par,'projectiveMaxMicroDt',Inf));
% 保证 M 个 micro steps 加 corrector steps 的物理时间不超过宏步。
dtm = min(dtm, dtMacro / max(M + corrector, 1));
if ~(isfinite(dtm) && dtm > 0)
    info.failed = true;
    info.failureReason = 'projective micro dt 非正或非有限';
    return;
end

Wprev = Wold; uprev = uold;
Wcur = Wold; ucur = uold;
ampInfoLast = struct();

for m = 1:M
    Wbefore = Wcur;
    ubefore = ucur;
    [rhoM, ~] = src.reconstruct_rho_1d(Wcur, ucur, par.eps, op, par.rhoReconstruction);
    HMraw = src.solve_H_fd_1d(D, Kx, rhoM, op);
    [HM, ~] = src.shift_effective_hamiltonian_1d(HMraw, par);
    [Wnext, unext, oneInfo] = src.wkb_one_step_frozenH_1d(Wcur, ucur, D, Kx, rhoM, HM, par.eps, dtm, op, par, variant);
    if oneInfo.failed
        info.failed = true;
        info.failureReason = ['projective micro step 失败: ' oneInfo.failureReason];
        info.ampInfo = oneInfo.ampInfo;
        return;
    end
    Wprev = Wbefore; uprev = ubefore;
    Wcur = Wnext; ucur = unext;
    ampInfoLast = oneInfo.ampInfo;
end

microTime = M * dtm;
rem = max(dtMacro - (M + corrector) * dtm, 0);
info.projectiveMicroStepsUsed = M;
info.projectiveMicroDt = dtm;
info.projectiveRemainingDt = rem;

% gauge-invariant 相位差分：去掉两个 micro 状态的常数 gauge，只外推形状。
up = uprev(:).' - max(uprev(:));
uc = ucur(:).' - max(ucur(:));
du = (uc - up) / dtm;
clipU = local_get_num(par,'projectiveDerivativeClipU',50);
du = min(max(du, -clipU), clipU);

unew = ucur + rem * du;
Wnew = Wcur;
info.projectiveModeCode = 1;

if any(strcmp(mode, {'u-logw','logw','full'}))
    % 可选：同时外推 log(W)。这更激进，可能更快，但也更容易在小 epsilon
    % 下制造指数跨度；因此不是默认。
    %
    % 实现细节：
    %   1. 只外推 log(W)，不直接外推 W，避免正性破坏；
    %   2. 对 d(log W)/dt 做限幅，防止宏步外推制造极端指数；
    %   3. 对外推后的 log(W) 再做绝对限幅，防止 exp 上溢；
    %   4. 若外推结果出现非有限值，则保留 micro-relaxed 的 Wcur 并标记失败。
    floorW = max(local_get_num(par,'WFloor',realmin), realmin);
    lwp = log(max(Wprev, floorW));
    lwc = log(max(Wcur, floorW));
    dlw = (lwc - lwp) / dtm;
    clipLW = local_get_num(par,'projectiveDerivativeClipLogW',50);
    dlw = min(max(dlw, -clipLW), clipLW);
    logWnew = lwc + rem * dlw;
    clipAbsLW = local_get_num(par,'projectiveLogWAbsClip',650);
    logWnew = min(max(logWnew, -clipAbsLW), clipAbsLW);
    if any(~isfinite(logWnew(:)))
        info.failed = true;
        info.failureReason = 'u-logw 外推得到非有限 log(W)';
        info.ampInfo = ampInfoLast;
        return;
    end
    Wnew = exp(logWnew);
    info.projectiveModeCode = 2;
end

% 常数 gauge 归一化，不改变重构密度。
[Wnew, unew] = src.normalize_gauge(Wnew, unew, par.eps, 'max');

% 外推后做少量 micro corrector，让 W 与新的 u/rho 重新靠近快流形。
for c = 1:corrector
    [rhoC, ~] = src.reconstruct_rho_1d(Wnew, unew, par.eps, op, par.rhoReconstruction);
    HCraw = src.solve_H_fd_1d(D, Kx, rhoC, op);
    [HC, ~] = src.shift_effective_hamiltonian_1d(HCraw, par);
    [Wnew, unew, oneInfo] = src.wkb_one_step_frozenH_1d(Wnew, unew, D, Kx, rhoC, HC, par.eps, min(dtm,dtMacro), op, par, variant);
    if oneInfo.failed
        info.failed = true;
        info.failureReason = ['projective corrector 失败: ' oneInfo.failureReason];
        info.ampInfo = oneInfo.ampInfo;
        return;
    end
    ampInfoLast = oneInfo.ampInfo;
end

if any(~isfinite(Wnew(:))) || any(~isfinite(unew(:)))
    info.failed = true;
    info.failureReason = 'projective 宏步结束后 W/u 出现 NaN/Inf';
end
info.ampInfo = ampInfoLast;
info.forceHistory = 1;
end


%==========================================================================
% 中间结果保存工具
%==========================================================================
function saveCtl = local_intermediate_save_init(par, tag)
%LOCAL_INTERMEDIATE_SAVE_INIT 初始化中间输出控制器。
%
%  1. full snapshot：直接保存到 outdir/snapshot_wkb_stepXXXXXXX_tT.mat，
%     适合保留若干中间时刻用于后处理。
%  2. latest checkpoint：保存到 outdir/checkpoint_wkb_latest.mat，每次覆盖，
%     适合中途停止后读取最近状态。
%  3. history checkpoint：保存到 outdir/history_wkb_latest.mat，只含 history
%     和摘要，文件较小，适合快速查看 residual 与 thetaWKB 曲线。
%
%  这些保存操作只记录当前 WKB 状态，不改变数值格式。

saveCtl = struct();
saveCtl.tag = tag;
saveCtl.outdir = utils.ensure_dir(par.outdir);
% 快照文件直接写入 outdir，不再创建 snapshots 子文件夹。
saveCtl.snapshotDir = saveCtl.outdir;

saveInitialSnapshot = local_get_bool(par,'saveInitialSnapshot',false);
saveCtl.snapshotEnabled = local_get_bool(par,'storeSnapshots',false) || ...
                          local_get_bool(par,'saveIntermediateSnapshots',false) || ...
                          saveInitialSnapshot;
if isfield(par,'saveIntermediate') && ~isempty(par.saveIntermediate)
    saveCtl.snapshotEnabled = local_get_bool(par,'saveIntermediate',saveCtl.snapshotEnabled) || saveInitialSnapshot;
end
saveCtl.checkpointEnabled = local_get_bool(par,'saveLatestCheckpoint',true);
saveCtl.historyEnabled = local_get_bool(par,'saveHistoryCheckpoint',true);

if isfield(par,'snapshotTimes') && ~isempty(par.snapshotTimes)
    saveCtl.snapshotTimes = sort(par.snapshotTimes(:).');
else
    saveCtl.snapshotTimes = [];
end
if saveInitialSnapshot
    saveCtl.snapshotTimes = unique([0, saveCtl.snapshotTimes]);
end
saveCtl.nextSnapshotIndex = 1;

saveCtl.snapshotEveryTime = local_get_num(par,'snapshotEveryTime',Inf);
saveCtl.snapshotEverySteps = round(local_get_num(par,'snapshotEverySteps',Inf));
if ~(isfinite(saveCtl.snapshotEveryTime) && saveCtl.snapshotEveryTime > 0)
    saveCtl.snapshotEveryTime = Inf;
end
if ~(isfinite(saveCtl.snapshotEverySteps) && saveCtl.snapshotEverySteps > 0)
    saveCtl.snapshotEverySteps = Inf;
end
saveCtl.nextPeriodicSnapshotTime = saveCtl.snapshotEveryTime;

saveCtl.checkpointEveryTime = local_get_num(par,'checkpointEveryTime',1.0);
saveCtl.checkpointEverySteps = round(local_get_num(par,'checkpointEverySteps',Inf));
if ~(isfinite(saveCtl.checkpointEveryTime) && saveCtl.checkpointEveryTime > 0)
    saveCtl.checkpointEveryTime = Inf;
end
if ~(isfinite(saveCtl.checkpointEverySteps) && saveCtl.checkpointEverySteps > 0)
    saveCtl.checkpointEverySteps = Inf;
end
saveCtl.nextCheckpointTime = saveCtl.checkpointEveryTime;

saveCtl.latestFile = local_solver_data_file(saveCtl.outdir, 'wkb', 'checkpoint', par, 'checkpointFile');
saveCtl.historyFile = local_solver_data_file(saveCtl.outdir, 'wkb', 'history', par, 'historyCheckpointFile');
saveCtl.lastSnapshotStep = -Inf;
saveCtl.lastCheckpointStep = -Inf;
saveCtl.savedCount = 0;
end

function [saveCtl, didSave] = local_maybe_save_intermediate(saveCtl, par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, stepInfo, massInfo)
%LOCAL_MAYBE_SAVE_INTERMEDIATE 按设定条件保存中间结果。
%  这个函数只在“已接受的宏步”之后调用，因此不会保存被拒绝的 trial step。

didSave = false;
tolT = 1e-12 + 1e-10*max(1, abs(t));

saveFull = false;
reason = '';
if saveCtl.snapshotEnabled
    if step == 1 && local_get_bool(par,'saveInitialAcceptedSnapshot',true)
        saveFull = true;
        reason = 'first';
    end
    while saveCtl.nextSnapshotIndex <= numel(saveCtl.snapshotTimes) && ...
            t + tolT >= saveCtl.snapshotTimes(saveCtl.nextSnapshotIndex)
        saveFull = true;
        reason = sprintf('time%.6g', saveCtl.snapshotTimes(saveCtl.nextSnapshotIndex));
        saveCtl.nextSnapshotIndex = saveCtl.nextSnapshotIndex + 1;
    end
    while isfinite(saveCtl.nextPeriodicSnapshotTime) && t + tolT >= saveCtl.nextPeriodicSnapshotTime
        saveFull = true;
        reason = sprintf('periodic%.6g', saveCtl.nextPeriodicSnapshotTime);
        saveCtl.nextPeriodicSnapshotTime = saveCtl.nextPeriodicSnapshotTime + saveCtl.snapshotEveryTime;
    end
    if isfinite(saveCtl.snapshotEverySteps) && saveCtl.snapshotEverySteps > 0 && step > 0 && ...
            mod(step, saveCtl.snapshotEverySteps) == 0 && step ~= saveCtl.lastSnapshotStep
        saveFull = true;
        reason = sprintf('step%d', step);
    end
end

if saveFull
    vars = local_make_save_vars(par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, stepInfo, massInfo);
    vars.saveReason = reason;
    fname = utils.snapshot_filename(saveCtl.outdir, 'wkb', step, t);
    local_save_struct_mat(fname, vars, par);
    saveCtl.lastSnapshotStep = step;
    saveCtl.savedCount = saveCtl.savedCount + 1;
    didSave = true;
end

saveLatest = false;
if saveCtl.checkpointEnabled
    if step == 1
        saveLatest = true;
    end
    while isfinite(saveCtl.nextCheckpointTime) && t + tolT >= saveCtl.nextCheckpointTime
        saveLatest = true;
        saveCtl.nextCheckpointTime = saveCtl.nextCheckpointTime + saveCtl.checkpointEveryTime;
    end
    if isfinite(saveCtl.checkpointEverySteps) && saveCtl.checkpointEverySteps > 0 && step > 0 && ...
            mod(step, saveCtl.checkpointEverySteps) == 0 && step ~= saveCtl.lastCheckpointStep
        saveLatest = true;
    end
end

if saveLatest
    vars = local_make_save_vars(par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, stepInfo, massInfo);
    vars.saveReason = 'latest-checkpoint';
    local_save_struct_mat(saveCtl.latestFile, vars, par);
    saveCtl.lastCheckpointStep = step;
    didSave = true;
end

if saveCtl.historyEnabled && (saveFull || saveLatest)
    histVars = struct();
    histVars.history = history;
    histVars.par = par;
    histVars.tag = tag;
    histVars.t = t;
    histVars.step = step;
    histVars.residual = residual;
    histVars.diagnostics = dg;
    histVars.dtInfo = dtInfo;
    histVars.resInfo = resInfo;
    if nargin >= 20 && isstruct(massInfo), histVars.massInfo = massInfo; end
    local_save_struct_mat(saveCtl.historyFile, histVars, par);
end
end

function [saveCtl, didSave] = local_force_save_latest(saveCtl, par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, stepInfo, massInfo)
%LOCAL_FORCE_SAVE_LATEST 强制写一次 latest checkpoint。
didSave = false;
if ~isstruct(saveCtl) || ~isfield(saveCtl,'latestFile')
    saveCtl = local_intermediate_save_init(par, tag);
end
if saveCtl.checkpointEnabled
    vars = local_make_save_vars(par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, stepInfo, massInfo);
    vars.saveReason = 'final-or-failed-checkpoint';
    local_save_struct_mat(saveCtl.latestFile, vars, par);
    didSave = true;
end
if saveCtl.historyEnabled
    histVars = struct('history',history,'par',par,'tag',tag,'t',t,'step',step, ...
        'residual',residual,'diagnostics',dg,'dtInfo',dtInfo,'resInfo',resInfo,'massInfo',massInfo);
    local_save_struct_mat(saveCtl.historyFile, histVars, par);
end
end

function vars = local_make_save_vars(par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, stepInfo, massInfo)
%LOCAL_MAKE_SAVE_VARS 组织保存变量。
vars = struct();
vars.W = W;
vars.u = u;
vars.n = n;
vars.rho = rho;
vars.H = H;
vars.D = D;
vars.Kx = Kx;
vars.grid = grid;
vars.op = op;
vars.par = par;
vars.tag = tag;
vars.t = t;
vars.step = step;
vars.residual = residual;
vars.history = history;
vars.diagnostics = dg;
vars.resInfo = resInfo;
vars.dtInfo = dtInfo;
vars.stepInfo = stepInfo;
if nargin >= 20 && isstruct(massInfo)
    vars.massInfo = massInfo;
end
vars.createdBy = 'run_wkb_timeap_1d intermediate output';
vars.createdAt = datestr(now, 31);
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
if isfield(par, 'massCorrectionFormula') && ~isempty(par.massCorrectionFormula)
    info.formula = char(par.massCorrectionFormula);
else
    info.formula = 'trapezoid';
end
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
    diag.mass_correction_scale = get_field(massInfo, 'scale', 1);
    diag.mass_trial = get_field(massInfo, 'trialMass', NaN);
    diag.mass_corrected = get_field(massInfo, 'correctedMass', massState.mass);
else
    diag.mass_correction_enabled = false;
    diag.mass_correction_applied = false;
    diag.mass_correction_scale = 1;
    diag.mass_trial = NaN;
    diag.mass_corrected = massState.mass;
end
end

function local_save_struct_mat(filename, vars, par)
%LOCAL_SAVE_STRUCT_MAT 保存 struct 到 mat 文件。
try
    if local_get_bool(par,'saveMatV73',false)
        save(filename, '-struct', 'vars', '-v7.3');
    else
        save(filename, '-struct', 'vars');
    end
catch ME
    warning('run_wkb_timeap_1d:IntermediateSaveFailed', ...
        '保存中间结果失败：%s。文件：%s', ME.message, filename);
end
end

%==========================================================================
% 工具函数
%==========================================================================
function par = local_preconfigure(par)
%LOCAL_PRECONFIGURE time-AP 路径的参数补齐。
%  timeIntegrator 本身已经是显式选择的 v13 实验补丁；但本函数仍然不
%  擅自改变 direct-log/laplace-hybrid 等用户设置，除非打开 smallEpsAutoPatches。

autoPatch = local_get_bool(par, 'smallEpsAutoPatches', false) || ...
            local_get_bool(par, 'enableV13AutoPatches', false);
if ~isfield(par,'rhoReconstruction') || isempty(par.rhoReconstruction)
    par.rhoReconstruction = 'direct-log';
end
if autoPatch && isfield(par,'eps') && par.eps <= 1e-4 && any(strcmpi(par.rhoReconstruction, {'direct','direct-log'}))
    par.rhoReconstruction = 'laplace-hybrid';
    if local_get_bool(par,'verbose',true)
        fprintf('[WKB-1D] smallEpsAutoPatches=true: 实验版时间积分器使用 laplace-hybrid 重构 rho。\n');
    end
end
if ~isfield(par,'reactionDiscretization') || isempty(par.reactionDiscretization)
    par.reactionDiscretization = 'implicit';
end
if autoPatch && (strcmpi(par.reactionDiscretization,'implicit') || strcmpi(par.reactionDiscretization,'patankar'))
    par.reactionDiscretization = 'logistic-exact';
end
if ~isfield(par,'ifMode') || isempty(par.ifMode)
    par.ifMode = 'none';
end
if autoPatch && strcmpi(par.ifMode,'none')
    par.ifMode = 'h-only';
end
if ~isfield(par,'hamiltonianGauge') || isempty(par.hamiltonianGauge)
    par.hamiltonianGauge = 'none';
end
if autoPatch && strcmpi(par.hamiltonianGauge,'none')
    par.hamiltonianGauge = 'min';
end
if ~isfield(par,'thetaRephase') || isempty(par.thetaRephase)
    par.thetaRephase = false;
end
end

function variant = local_normalize_variant(name)
%LOCAL_NORMALIZE_VARIANT 统一振幅格式名称；保持 density-compatible 的旧语义。
name = lower(char(name));
switch name
    case {'wkb-if-split','wkb-if','balanced-wkb-if'}
        variant = 'wkb-if-split';
    case {'split-if','balanced-split','if-split'}
        variant = 'split-if';
    case {'split','split-wkb'}
        variant = 'split';
    case {'density-compatible','dc','density-compatible-semidiscrete','dc-semidiscrete','density-compatible-old','dc-old'}
        variant = 'density-compatible-semidiscrete';
    case {'density-compatible-if','dc-if','density-if','integrating-factor'}
        variant = 'wkb-if-split';
    otherwise
        error('Unknown amplitudeVariant "%s".', name);
end
end

function n = local_every_steps(par, fieldSteps, fieldTime, defaultTime)
if isfield(par, fieldSteps) && ~isempty(par.(fieldSteps)) && isnumeric(par.(fieldSteps)) && isscalar(par.(fieldSteps)) && par.(fieldSteps) > 0
    n = max(1, round(par.(fieldSteps)));
elseif isfield(par, fieldTime) && ~isempty(par.(fieldTime)) && isnumeric(par.(fieldTime)) && isscalar(par.(fieldTime)) && par.(fieldTime) > 0
    n = max(1, round(par.(fieldTime) / par.dt));
else
    n = max(1, round(defaultTime / par.dt));
end
end

function state = local_make_plot_state(W, u, n, rho, H, Kx, op, par, t, step, residual, dg)
state = struct();
state.W = real(W);
state.u = real(u);
state.n = real(n);
state.rho = real(rho);
state.H = real(H);
state.Kx = real(Kx);
state.x = op.x;
state.theta = op.theta;
state.P = real(dg.P);
state.theta_peak = dg.theta_wkb;
state.theta_m = par.theta_m;
state.sigma_theta = dg.sigma_theta;
state.gammaR = dg.gammaR;
state.rho_max = dg.rho_max;
state.residual = residual;
state.t = t;
state.T = par.T;
state.step = step;
state.saveFrame = false;
end

function info = local_progress_info(dg, resInfo, dtInfo, timeIntegrator)
if nargin < 2 || isempty(resInfo), resInfo = struct(); end
if nargin < 3 || isempty(dtInfo), dtInfo = struct(); end
if nargin < 4, timeIntegrator = ''; end
info = struct();
if isstruct(dg)
    info.theta_wkb = get_field(dg,'theta_wkb',NaN);
    info.theta_direct = get_field(dg,'theta_direct',NaN);
    info.err_wkb_to_m = get_field(dg,'err_wkb_to_m',NaN);
    info.sigma_theta = get_field(dg,'sigma_theta',NaN);
    info.gammaR = get_field(dg,'gammaR',NaN);
    info.rho_max = get_field(dg,'rho_max',NaN);
    info.Hmin = get_field(dg,'Hmin',NaN);
    info.res_phase_steady = get_field(dg,'res_phase_steady',NaN);
    info.res_density_steady = get_field(dg,'res_density_steady',NaN);
    info.res_amp_w_steady = get_field(dg,'res_amplitude_steady',get_field(dg,'res_amp_w_steady',NaN));
    info.phaseCFL = get_field(dg,'phaseCFL',NaN);
    info.jumpEps = get_field(dg,'neighborJumpOverEps',NaN);
    info.peakGrid = get_field(dg,'peakWidthOverGrid',NaN);
    info.logQMax = get_field(dg,'logQMax',NaN);
    info.seedLogRange = get_field(dg,'seedLogRange',NaN);
    info.reactionScaleMax = get_field(dg,'reactionScaleMax',NaN);
    info.implicitIter = get_field(dg,'implicitIter',NaN); %#ok<STRNU>
    info.implicitRelChange = get_field(dg,'implicitRelChange',NaN); %#ok<STRNU>
    info.projectiveMicroStepsUsed = get_field(dg,'projectiveMicroStepsUsed',NaN); %#ok<STRNU>
    info.projectiveMicroDt = get_field(dg,'projectiveMicroDt',NaN); %#ok<STRNU>
end
if isstruct(resInfo)
    if isfield(resInfo,'res_n'), info.res_n = resInfo.res_n; end
    if isfield(resInfo,'res_u'), info.res_u = resInfo.res_u; end
    if isfield(resInfo,'res_W_legacy'), info.res_W_legacy = resInfo.res_W_legacy; end
end
if isstruct(dtInfo)
    info.dt = get_field(dtInfo,'dt',NaN);
    info.lambdaR = get_field(dtInfo,'lambdaR',NaN);
    info.lambdaH = get_field(dtInfo,'lambdaH',NaN);
    info.ifLogRangeEstimate = get_field(dtInfo,'ifLogRangeEstimate',NaN);
    if isfield(dtInfo,'numRetries'), info.numRetries = dtInfo.numRetries; end
end
if ~isempty(timeIntegrator)
    info.note = ['timeIntegrator=' char(timeIntegrator)];
end
end

function c = local_integrator_code(s)
s = lower(char(s));
if any(strcmp(s, {'implicit-coupled','fully-implicit-coupled','coupled-implicit'}))
    c = 1;
elseif any(strcmp(s, {'projective','micro-macro','projective-micro'}))
    c = 3;
else
    c = 0;
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

function v = get_field(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
end

function v = local_get_num(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
end


function filename = local_solver_data_file(outdir, solverName, kind, par, fieldName)
filename = utils.solver_file(outdir, solverName, kind);
if isstruct(par) && isfield(par, fieldName) && ~isempty(par.(fieldName))
    customName = char(par.(fieldName));
    if utils.is_absolute_path(customName)
        filename = customName;
    else
        filename = fullfile(outdir, customName);
    end
end
end

function v = local_get_string(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    v = char(s.(name));
end
end

function tf = local_get_bool(s, name, defaultValue)
tf = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end
