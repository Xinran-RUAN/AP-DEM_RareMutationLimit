function result = run_wkb_1d(par)
%RUN_WKB_1D  运行一维有限差分 WKB 求解器。
%
%  主未知量始终是 WKB 变量 (W,u)：
%     rho = Reconstruction[W,u],   H = principal eigenvalue,
%     u^{n+1} from HJ equation,    W^{n+1} from WKB amplitude equation.
%
%  本版本针对小 epsilon 增加了三类保护：
%   1. phase-peak gauge：每步把相位的二次插值峰值规范到 0，防止
%      Laplace 重构里的 exp(u*/eps) 上溢；
%   2. WKB-IF 的 RHS log 监测：若 qW 的 log 尺度过大，拒绝本步并缩小 dt；
%   3. 详细中文诊断：记录相位残差、振幅残差、rho 重构、q 指数、峰宽/网格比。
%
%  这些处理只是在 WKB 变量内部做 gauge/时间步/重构稳健化，不是直接求解 n 方程。

% v11: 可选的时间 AP 实验版求解器。
%   当 par.timeIntegrator 取 'implicit-coupled' 或 'projective' 时，
%   直接转入专门的实验版主循环。这样保留原来的 frozen-H 时间推进路径，
%   同时可以独立测试“全隐式耦合 H^{n+1}”和“projective/micro-macro”两种方案。
if nargin >= 1 && isstruct(par) && isfield(par,'timeIntegrator') && ~isempty(par.timeIntegrator)
    ti__ = lower(char(par.timeIntegrator));
    if any(strcmp(ti__, {'implicit-coupled','fully-implicit-coupled','coupled-implicit', ...
                         'projective','micro-macro','projective-micro'}))
        result = src.run_wkb_timeap_1d(par);
        return;
    end
end

par = local_preconfigure_small_eps(par);
[W, u, D, Kx, grid] = model.initial_wkb_1d(par);
op = src.build_operators_1d(grid);

% 小 epsilon 下的关键稳健化：theta 方向再分解。
% 它精确保持 n=W*exp(u/eps) 不变，但把 W 中积累的 theta 指数尺度
% 吸收到相位 u 中，避免后续 WKB-IF 右端 qW 出现巨大 log 范围。
[W, u, initRephaseInfo] = src.rebalance_wkb_theta_gauge_1d(W, u, par.eps, op, par);

par.outdir = utils.ensure_dir(par.outdir);
tag = utils.build_tag(par);
prog = utils.progress_init('WKB-1D', par);
livePlotEverySteps = local_liveplot_every_steps(par);
lp = utils.liveplot_init_1d(par, 'wkb');

% 中间结果保存控制器。普通 frozen-H 路径也保存 snapshots/checkpoint，
% 与 time-AP/projective 路径保持一致。
saveCtl = utils.intermediate_save_init(par, tag, 'wkb');
history = [];
historyEverySteps = local_history_every_steps(par);
lastInfo = struct();
lastPlotState = struct();

variant = local_normalize_variant(par.amplitudeVariant);

t = 0;
step = 0;
residual = Inf;
failed = false;
failureReason = '';
resInfo = struct();
dtInfo = struct();

% 质量平衡修正状态。massState 始终对应最近一个已接受状态。WKB 中质量修正只缩放 W，u 不变。
[rhoInitialForMass, ~] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
massState = src.mass_balance_stats_1d(rhoInitialForMass, Kx, op);
massInfoCurrent = local_make_initial_mass_info(massState, par);

% 可选保存初始时刻 t=0。旧版若 snapshotTimes 含 0 或 saveInitialSnapshot=true，
% 用户通常希望后处理能看到初值；v13 原实现只在 step=1 后保存，容易丢失 t=0。
[initialSaveDue, ~, ~] = utils.intermediate_save_due(saveCtl, t, step);
if initialSaveDue
    try
        [rho0, n0] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
        H0raw = src.solve_H_fd_1d(D, Kx, rho0, op);
        [H0, hGaugeInfo0] = src.shift_effective_hamiltonian_1d(H0raw, par);
        dg0 = src.compute_diagnostics_1d(W, u, D, Kx, rho0, H0, par.eps, op, par, variant);
        dg0 = local_merge_structs(dg0, hGaugeInfo0, initRephaseInfo);
        dg0 = local_attach_mass_info(dg0, massInfoCurrent, massState);
        vars0 = local_make_save_vars(par, tag, t, step, W, u, n0, rho0, H0, D, Kx, grid, op, residual, history, dg0, resInfo, dtInfo, variant, massInfoCurrent);
        [saveCtl, didSave0] = utils.intermediate_save_step(saveCtl, par, vars0);
        if didSave0 && par.verbose && isfield(par,'printIntermediateSaveMessage') && par.printIntermediateSaveMessage
            fprintf('[WKB-1D] 已保存初始中间结果: step=0, t=0\n');
        end
    catch ME
        warning('run_wkb_1d:InitialSnapshotFailed', '初始 WKB 快照保存失败：%s', ME.message);
    end
end

while t < par.T && step < par.maxSteps
    retry = 0;
    maxRetry = local_get_field(par, 'maxStepRetries', 12);
    retryFactor = local_get_field(par, 'retryDtFactor', 0.5);
    dtUpper = Inf;
    accepted = false;
    rejectReason = '';

    while ~accepted
        % 每次重试都从同一个已接受状态 (W,u,t) 出发；不会污染当前解。
        [rho, ~, rhoInfo0] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
        if any(~isfinite(rho(:)))
            failed = true;
            failureReason = '当前状态的 rho 重构已经出现 NaN/Inf';
            break;
        end
        Hraw = src.solve_H_fd_1d(D, Kx, rho, op);
        if any(~isfinite(Hraw(:)))
            failed = true;
            failureReason = '当前状态的 H 特征值计算已经出现 NaN/Inf';
            break;
        end
        % 对 H 做标量 gauge 平移。后续 phase 与 amplitude 必须使用同一个 H。
        % 默认 min(Huse)=0，可去掉 H 的公共偏移，避免 exp(dt*H/eps) 产生无意义刚性。
        [H, hGaugeInfo] = src.shift_effective_hamiltonian_1d(Hraw, par);

        [dtStep, dtInfo] = src.select_wkb_time_step_1d(par, t, W, u, D, Kx, rho, H, op, variant);
        dtStep = min(dtStep, dtUpper);
        % select_wkb_time_step_1d 的 dtInfo 是按原候选步长计算的；若这里因拒步
        % dtUpper 再次缩小，需要同步更新与 dt 成正比的诊断量。
        dtInfo.dt = dtStep;
        dtInfo.dtOverEps = dtStep / par.eps;
        if isfield(dtInfo,'rplusMax'), dtInfo.lambdaR = dtStep/par.eps * dtInfo.rplusMax; end
        if isfield(dtInfo,'HabsMax'),  dtInfo.lambdaH = dtStep/par.eps * dtInfo.HabsMax; end
        if isfield(dtInfo,'phaseSpeedMax'), dtInfo.phaseCFL = dtInfo.phaseSpeedMax * dtStep / op.dtheta; end
        dtInfo.ifLogRangeEstimate = dtStep/par.eps * max(max(H(:))-min(H(:)),0);
        if dtStep <= 0 || ~isfinite(dtStep)
            failed = true;
            failureReason = 'dtStep 非正或非有限';
            break;
        end

        uOld = u;
        WOld = W;
        ampInfo = struct();
        gaugeInfo = struct();
        postGaugeInfo = struct();
        rephaseInfo = struct();
        massInfoTrial = massInfoCurrent;
        massStateTrial = massState;

        %-------------------------
        % 相位方程。相位方程本身只依赖 u 的导数，因此后续常数 gauge 不影响 HJ 结构。
        %-------------------------
        uRaw = src.step_phase_1d(uOld, H, par.eps, dtStep, op, par.phaseHamiltonian, par.lfAlpha);

        %-------------------------
        % 振幅方程。所有分支仍然更新 W，而不是求解 n 方程。
        %-------------------------
        switch variant
            case 'wkb-if-split'
                [WTrial, uTrial, ampInfo] = src.update_w_wkb_if_split_fd_1d(WOld, uOld, uRaw, D, Kx, rho, H, par.eps, dtStep, op, par);
                gaugeInfo = struct('phaseGaugeShift', local_get_field(ampInfo, 'ifGaugeShift', 0));
                postGaugeInfo = struct();
            case 'split-if'
                [uAmp, ~, gaugeInfo] = src.choose_phase_gauge_1d(uOld, uRaw, par.eps, par);
                [WStep, ampInfo] = src.update_w_split_if_fd_1d(WOld, uOld, uAmp, D, Kx, rho, H, par.eps, dtStep, op, par);
                [WTrial, uTrial, postGaugeInfo] = src.stabilize_wkb_gauge_1d(WStep, uAmp, par.eps, par);
            case 'split'
                [WStep, ampInfo] = src.update_w_split_fd_1d(WOld, uRaw, D, Kx, rho, H, par.eps, dtStep, op, par);
                [WTrial, uTrial] = src.normalize_gauge(WStep, uRaw, par.eps, 'max');
            case 'density-compatible-semidiscrete'
                [WStep, ampInfo] = src.update_w_density_compatible_fd_1d(WOld, uRaw, D, Kx, rho, H, par.eps, dtStep, op, par);
                [WTrial, uTrial] = src.normalize_gauge(WStep, uRaw, par.eps, 'max');
            case 'density-compatible-if'
                [WStep, ampInfo] = src.update_w_density_if_fd_1d(WOld, uOld, uRaw, D, Kx, rho, par.eps, dtStep, op, par);
                [WTrial, uTrial] = src.normalize_gauge(WStep, uRaw, par.eps, 'max');
        end

        % v10 默认不再作 theta 依赖 rephase。
        % 原因：theta 依赖 rephase 会把 log(W) 吸收到 u 中，虽然保持 n 不变，
        % 但会改变 u_theta，从而可能把相位峰人为锁在初始 trait 附近。
        rephaseInfo = struct('thetaRephaseApplied',0);

        % 先检查振幅更新是否明确要求重试。
        rejectThis = false;
        if isstruct(ampInfo) && isfield(ampInfo,'needsRetry') && ampInfo.needsRetry
            rejectThis = true;
            rejectReason = ampInfo.retryReason;
        end
        if any(~isfinite(WTrial(:))) || any(~isfinite(uTrial(:)))
            rejectThis = true;
            rejectReason = 'WTrial/uTrial 出现 NaN 或 Inf';
        end

        % 每一步可选进行质量平衡投影。WKB 情形下只缩放 W，保持 u 和峰值位置不变。
        if ~rejectThis && local_mass_correction_requested(par)
            try
                [rhoTrialForMass, ~] = src.reconstruct_rho_1d(WTrial, uTrial, par.eps, op, par.rhoReconstruction);
                [massScale, massInfoTrial] = src.mass_balance_correction_1d( ...
                    rhoTrialForMass, Kx, op, par.eps, dtStep, ...
                    massState.mass, massState.reactionIntegral, par);
                if massInfoTrial.applied
                    WTrial = massScale * WTrial;
                    rhoTrialForMass = massScale * rhoTrialForMass;
                end
                massStateTrial = src.mass_balance_stats_1d(rhoTrialForMass, Kx, op);
                massInfoTrial = local_finalize_mass_info(massInfoTrial, massStateTrial);
            catch ME
                rejectThis = true;
                rejectReason = ['质量平衡修正失败: ' ME.message];
            end
        end

        % 计算 trial residual。若 residual 或 rho/H 诊断已经非有限，也拒绝本步。
        residualMode = lower(local_get_string(par, 'residualMode', 'wkb-steady'));
        if any(strcmp(residualMode, {'wkb-steady','steady','steady-wkb'}))
            updateModeForResidual = 'scaled-gauge-invariant';
        else
            updateModeForResidual = residualMode;
        end
        [updateResidual, resInfoTrial] = src.compute_wkb_residual_1d(WTrial, uTrial, WOld, uOld, par.eps, dtStep, op, updateModeForResidual);

        steadyInfoTrial = struct();
        if any(strcmp(residualMode, {'wkb-steady','steady','steady-wkb'}))
            [rhoRes, ~, rhoInfoTrial] = src.reconstruct_rho_1d(WTrial, uTrial, par.eps, op, par.rhoReconstruction);
            rhoLogMaxTrial = local_get_field(rhoInfoTrial,'logrhoMax',NaN);
            rhoCapCountTrial = local_get_field(rhoInfoTrial,'numCappedHigh',local_get_field(rhoInfoTrial,'numCappedRhoHigh',0));
            rhoLogRetryMax = local_get_field(par,'rhoLogRetryMax',120);
            if any(~isfinite(rhoRes(:)))
                rejectThis = true;
                rejectReason = 'trial rho 重构出现 NaN 或 Inf';
            elseif (isfinite(rhoLogMaxTrial) && rhoLogMaxTrial > rhoLogRetryMax) || rhoCapCountTrial > 0
                rejectThis = true;
                rejectReason = sprintf('trial rho 的 log 尺度过大: logRhoMax=%.3g, capped=%g', rhoLogMaxTrial, rhoCapCountTrial);
            else
                HRes = src.solve_H_fd_1d(D, Kx, rhoRes, op);
                if any(~isfinite(HRes(:)))
                    rejectThis = true;
                    rejectReason = 'trial H 计算出现 NaN 或 Inf';
                else
                    steadyInfoTrial = src.compute_wkb_steady_residual_1d(WTrial, uTrial, D, Kx, rhoRes, HRes, par.eps, op, par, variant);
                    updateResidual = steadyInfoTrial.res_wkb_steady;
                    steadyInfoTrial.rhoReconstructLogMax = local_get_field(rhoInfoTrial,'logrhoMax',NaN);
                    steadyInfoTrial.rhoReconstructNumCappedHigh = local_get_field(rhoInfoTrial,'numCappedHigh',local_get_field(rhoInfoTrial,'numCappedRhoHigh',0));
                end
            end
        end

        if ~isfinite(updateResidual)
            rejectThis = true;
            if isempty(rejectReason), rejectReason = 'trial residual 非有限'; end
        end

        % 用户可要求“失败即缩步重试”。默认开启，尤其用于 eps<<1 的排查。
        useRetry = local_get_bool(par, 'enableStepRetry', false);
        dtMinRetry = local_get_field(par, 'dtMinRetry', max(1e-12, 1e-8*par.dt));
        if rejectThis && useRetry && retry < maxRetry && dtStep*retryFactor >= dtMinRetry
            retry = retry + 1;
            dtUpper = dtStep * retryFactor;
            if par.verbose && (retry == 1 || retry == maxRetry)
                fprintf('[WKB-1D] 本步被拒绝，缩小 dt 重试: t=%.6g, old dt=%.3e, new dt<=%.3e, 原因: %s\n', ...
                    t, dtStep, dtUpper, rejectReason);
            end
            continue;
        end

        % 无法再重试时，若仍失败则停止；否则接受本步。
        if rejectThis
            failed = true;
            failureReason = rejectReason;
            break;
        end

        if ~local_mass_correction_requested(par)
            [rhoTrialNoCorr, ~] = src.reconstruct_rho_1d(WTrial, uTrial, par.eps, op, par.rhoReconstruction);
            massStateTrial = src.mass_balance_stats_1d(rhoTrialNoCorr, Kx, op);
            massInfoTrial = local_finalize_mass_info(massInfoTrial, massStateTrial);
        end

        accepted = true;
        W = WTrial;
        u = uTrial;
        residual = updateResidual;
        resInfo = resInfoTrial;
        steadyInfo = steadyInfoTrial;
        massState = massStateTrial;
        massInfoCurrent = massInfoTrial;
        dtInfo.dt = dtStep;
        dtInfo.dtLimitedByRetry = retry > 0;
        dtInfo.numRetries = retry;
        dtInfo.rhoReconstructLogMaxOld = local_get_field(rhoInfo0,'logrhoMax',NaN);
    end

    if failed
        residual = Inf;
        if par.verbose
            fprintf('\n[WKB-1D] 警告：检测到数值失败，提前停止。原因: %s\n', failureReason);
        end
        break;
    end

    t = t + dtStep;
    step = step + 1;

    % 小 epsilon 监测量：每步便宜计算，history 或保存时写入文件。
    monInfo = local_small_eps_monitor(W, u, uOld, dtStep, op, par, ampInfo, gaugeInfo, postGaugeInfo);

    historyDue = (step == 1 || mod(step, historyEverySteps) == 0 || residual <= par.tol || retry > 0);
    plotDue = lp.enabled && (step == 1 || mod(step, livePlotEverySteps) == 0 || residual <= par.tol);
    [saveDue, ~, ~] = utils.intermediate_save_due(saveCtl, t, step);

    if historyDue || plotDue || saveDue
        [rhoNow, nNow] = src.reconstruct_rho_1d(W, u, par.eps, op, par.rhoReconstruction);
        HNowRaw = src.solve_H_fd_1d(D, Kx, rhoNow, op);
        [HNow, hGaugeInfoNow] = src.shift_effective_hamiltonian_1d(HNowRaw, par);
        massState = src.mass_balance_stats_1d(rhoNow, Kx, op);
        massInfoCurrent = local_finalize_mass_info(massInfoCurrent, massState);
        dg = src.compute_diagnostics_1d(W, u, D, Kx, rhoNow, HNow, par.eps, op, par, variant);
        dg = local_merge_structs(dg, hGaugeInfoNow, steadyInfo, ampInfo, gaugeInfo, postGaugeInfo, rephaseInfo, monInfo);
        dg = local_attach_mass_info(dg, massInfoCurrent, massState);
        dg.numRetries = retry;

        if historyDue
            history = utils.append_history(history, t, residual, dg, resInfo, dtInfo);
        end
        lastInfo = local_progress_info(dg, resInfo, dtInfo);
        lastPlotState = local_make_plot_state(W, u, nNow, rhoNow, HNow, Kx, op, par, t, step, residual, dg);

        if saveDue
            vars = local_make_save_vars(par, tag, t, step, W, u, nNow, rhoNow, HNow, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, variant, massInfoCurrent);
            [saveCtl, didSaveIntermediate] = utils.intermediate_save_step(saveCtl, par, vars);
            if didSaveIntermediate && par.verbose && isfield(par,'printIntermediateSaveMessage') && par.printIntermediateSaveMessage
                fprintf('[WKB-1D] 已保存中间结果: step=%d, t=%.6g\n', step, t);
            end
        end
    end

    prog = utils.progress_update(prog, step, t, residual, lastInfo);

    if plotDue && ~isempty(fieldnames(lastPlotState))
        lastPlotState.saveFrame = isfield(par, 'livePlotSave') && par.livePlotSave;
        lp = utils.liveplot_update_1d(lp, 'wkb', lastPlotState);
        if isfield(par, 'livePlotPause') && par.livePlotPause > 0
            pause(par.livePlotPause);
        end
    end

    if par.stopByResidual && isfinite(residual) && residual <= par.tol
        break;
    end
end

% 最终保存。若最终 rho/H 仍失败，也把失败原因写入 result，而不是误报收敛。
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
    diagInfo = struct('theta_wkb',NaN,'err_wkb_to_m',NaN,'sigma_theta',NaN,'rho_max',Inf,'gammaR',Inf, ...
        'failureReason', ['final diagnostics failed: ' ME.message]);
    residual = Inf;
end
if failed
    diagInfo.failed = true;
    diagInfo.failureReason = failureReason;
else
    diagInfo.failed = false;
    diagInfo.failureReason = '';
end

result = struct('W',W,'u',u,'n',n,'rho',rho,'H',H,'D',D,'Kx',Kx,'grid',grid, ...
                'op',op,'par',par,'t',t,'step',step,'residual',residual, ...
                'history',history,'diagnostics',diagInfo,'tag',tag,'variant',variant, ...
                'massInfo',massInfoCurrent);

% 无论正常结束、达到终止条件、失败退出，最后都补写一次 latest checkpoint。
try
    finalVars = local_make_save_vars(par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, diagInfo, struct(), struct(), variant, massInfoCurrent);
    [saveCtl, ~] = utils.intermediate_save_force_latest(saveCtl, par, finalVars, 'final-or-failed-checkpoint'); %#ok<NASGU>
catch ME
    warning('run_wkb_1d:FinalCheckpointFailed', '最终 WKB checkpoint 保存失败：%s', ME.message);
end

utils.progress_finish(prog, step, t, residual, local_progress_info(diagInfo, struct(), struct()));
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




function vars = local_make_save_vars(par, tag, t, step, W, u, n, rho, H, D, Kx, grid, op, residual, history, dg, resInfo, dtInfo, variant, massInfo)
% 组织普通 WKB/frozen-H 路径的完整中间输出变量。
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
vars.variant = variant;
if nargin >= 20 && isstruct(massInfo)
    vars.massInfo = massInfo;
end
vars.createdBy = 'run_wkb_1d intermediate output';
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
if nargin < 3 || ~isstruct(massState)
    massState = struct('mass',NaN,'reactionIntegral',NaN);
end
diag.mass_total = massState.mass;
diag.mass_balance_rhs = massState.reactionIntegral;
if isstruct(massInfo)
    diag.mass_correction_enabled = local_get_bool(massInfo, 'enabled', false);
    diag.mass_correction_applied = local_get_bool(massInfo, 'applied', false);
    diag.mass_correction_scale = local_get_field(massInfo, 'scale', 1);
    diag.mass_trial = local_get_field(massInfo, 'trialMass', NaN);
    diag.mass_corrected = local_get_field(massInfo, 'correctedMass', massState.mass);
else
    diag.mass_correction_enabled = false;
    diag.mass_correction_applied = false;
    diag.mass_correction_scale = 1;
    diag.mass_trial = NaN;
    diag.mass_corrected = massState.mass;
end
end

function par = local_preconfigure_small_eps(par)
%LOCAL_PRECONFIGURE_SMALL_EPS 只补齐缺失字段；是否自动改写由开关控制。
%
%  旧版兼容原则：默认不因为 eps 很小就擅自改变 rhoReconstruction、
%  reactionDiscretization 或 Hamiltonian gauge。若用户显式打开
%      par.smallEpsAutoPatches = true
%  或
%      par.enableV13AutoPatches = true,
%  才应用 v13 推荐的小 epsilon 稳健默认值。

autoPatch = local_get_bool(par, 'smallEpsAutoPatches', false) || ...
            local_get_bool(par, 'enableV13AutoPatches', false);

if ~isfield(par,'rhoReconstruction') || isempty(par.rhoReconstruction)
    par.rhoReconstruction = 'direct-log';
end
if ~isfield(par,'reactionDiscretization') || isempty(par.reactionDiscretization)
    par.reactionDiscretization = 'implicit';
end
if ~isfield(par,'ifMode') || isempty(par.ifMode)
    par.ifMode = 'none';
end
if ~isfield(par,'hamiltonianGauge') || isempty(par.hamiltonianGauge)
    par.hamiltonianGauge = 'none';
end
if ~isfield(par,'thetaRephase') || isempty(par.thetaRephase)
    par.thetaRephase = false;
end
if ~isfield(par,'thetaRephaseStatistic') || isempty(par.thetaRephaseStatistic)
    par.thetaRephaseStatistic = 'median';
end
if ~isfield(par,'thetaRephaseSmoothPasses') || isempty(par.thetaRephaseSmoothPasses)
    par.thetaRephaseSmoothPasses = 0;
end

if ~autoPatch
    return;
end

% v13 推荐设置只在显式打开 autoPatch 后生效。
threshold = 1e-4;
mode = lower(char(par.rhoReconstruction));
if isfield(par,'eps') && par.eps <= threshold && any(strcmp(mode, {'direct','direct-log'}))
    par.rhoReconstruction = 'laplace-hybrid';
    if local_get_bool(par,'verbose',true)
        fprintf('[WKB-1D] smallEpsAutoPatches=true: rhoReconstruction 从 direct-log 切换为 laplace-hybrid。\n');
    end
end
if strcmpi(par.reactionDiscretization,'patankar') || strcmpi(par.reactionDiscretization,'implicit')
    par.reactionDiscretization = 'logistic-exact';
end
if strcmpi(par.ifMode,'none')
    par.ifMode = 'h-only';
end
if strcmpi(par.hamiltonianGauge,'none')
    par.hamiltonianGauge = 'min';
end
end

function failInfo = local_failure_check(W, u, residual, par)
failInfo = struct('failed', false, 'failureReason', '');
if any(~isfinite(W(:))) || any(~isfinite(u(:))) || ~isfinite(residual)
    failInfo.failed = true;
    failInfo.failureReason = 'W/u/residual 出现 NaN 或 Inf';
    return;
end
thrW = local_get_field(par, 'WBlowupThreshold', 1e200);
if max(abs(W(:))) > thrW
    failInfo.failed = true;
    failInfo.failureReason = 'W 超过 WBlowupThreshold';
    return;
end
thru = local_get_field(par, 'uBlowupThreshold', 1e100);
if max(abs(u(:))) > thru
    failInfo.failed = true;
    failInfo.failureReason = 'u 超过 uBlowupThreshold';
    return;
end
end

function variant = local_normalize_variant(name)
%LOCAL_NORMALIZE_VARIANT 统一振幅格式名称。
%  为了旧版兼容，'density-compatible'/'dc' 明确指向旧半离散矩阵版本；
%  只有显式写 'density-compatible-if'/'dc-if' 或 'wkb-if-split' 时才使用 v13 IF 格式。
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

function historyEverySteps = local_history_every_steps(par)
if isfield(par, 'historyEverySteps') && ~isempty(par.historyEverySteps) && ...
        isnumeric(par.historyEverySteps) && isscalar(par.historyEverySteps) && par.historyEverySteps > 0
    historyEverySteps = max(1, round(par.historyEverySteps));
elseif isfield(par, 'historyEveryTime') && ~isempty(par.historyEveryTime) && ...
        isnumeric(par.historyEveryTime) && isscalar(par.historyEveryTime) && par.historyEveryTime > 0
    historyEverySteps = max(1, round(par.historyEveryTime / par.dt));
else
    historyEverySteps = max(1, round(0.1 / par.dt));
end
end

function livePlotEverySteps = local_liveplot_every_steps(par)
if isfield(par, 'livePlotEverySteps') && ~isempty(par.livePlotEverySteps) && ...
        isnumeric(par.livePlotEverySteps) && isscalar(par.livePlotEverySteps) && par.livePlotEverySteps > 0
    livePlotEverySteps = max(1, round(par.livePlotEverySteps));
elseif isfield(par, 'livePlotEveryTime') && ~isempty(par.livePlotEveryTime) && ...
        isnumeric(par.livePlotEveryTime) && isscalar(par.livePlotEveryTime) && par.livePlotEveryTime > 0
    livePlotEverySteps = max(1, round(par.livePlotEveryTime / par.dt));
else
    livePlotEverySteps = max(1, round(1 / par.dt));
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

function info = local_progress_info(dg, resInfo, dtInfo)
if nargin < 2 || isempty(resInfo), resInfo = struct(); end
if nargin < 3 || isempty(dtInfo), dtInfo = struct(); end
info = struct();
if ~isempty(dg)
    info.theta_wkb = local_get_field(dg,'theta_wkb',NaN);
    info.theta_direct = local_get_field(dg,'theta_direct',NaN);
    info.err_wkb_to_m = local_get_field(dg,'err_wkb_to_m',NaN);
    info.sigma_theta = local_get_field(dg,'sigma_theta',NaN);
    info.gammaR = local_get_field(dg,'gammaR',NaN);
    info.rho_max = local_get_field(dg,'rho_max',NaN);
    info.mass_total = local_get_field(dg,'mass_total',NaN);
    info.mass_correction_scale = local_get_field(dg,'mass_correction_scale',NaN);
    info.Hmin = local_get_field(dg,'Hmin',NaN);
    info.qLogAbsMax = local_get_field(dg,'qLogAbsMax',NaN);
    info.logQMax = local_get_field(dg,'logQMax', info.qLogAbsMax);
    info.phaseCFL = local_get_field(dg,'phaseCFL',NaN);
    info.jumpEps = local_get_field(dg,'neighborJumpOverEps',NaN);
    info.peakGrid = local_get_field(dg,'peakWidthOverGrid',NaN);
    info.res_phase_steady = local_get_field(dg,'res_phase_steady',NaN);
    info.res_density_steady = local_get_field(dg,'res_density_steady',NaN);
    info.res_amplitude_steady = local_get_field(dg,'res_amplitude_steady',local_get_field(dg,'res_amp_w_steady',NaN));
    info.res_amp_w_steady = info.res_amplitude_steady;
    info.ifGaugeShift = local_get_field(dg,'ifGaugeShift',NaN);
    info.numRetries = local_get_field(dg,'numRetries',NaN);
    info.seedLogRange = local_get_field(dg,'seedLogRange',NaN);
    info.seedLogMax = local_get_field(dg,'seedLogMax',NaN);
    info.reactionScaleMax = local_get_field(dg,'reactionScaleMax',NaN);
    info.rhoLogMax = local_get_field(dg,'rhoReconstructLogMax',NaN);
    info.logWRangeAfter = local_get_field(dg,'logWRangeAfter',NaN);
    info.thetaRephaseARange = local_get_field(dg,'thetaRephaseARange',NaN);
end
if isstruct(resInfo)
    if isfield(resInfo,'res_n'), info.res_n = resInfo.res_n; end
    if isfield(resInfo,'res_u'), info.res_u = resInfo.res_u; end
    if isfield(resInfo,'res_W_legacy'), info.res_W_legacy = resInfo.res_W_legacy; end
end
if isstruct(dtInfo)
    if isfield(dtInfo,'dt'), info.dt = dtInfo.dt; end
    if isfield(dtInfo,'lambdaR'), info.lambdaR = dtInfo.lambdaR; end
    if isfield(dtInfo,'lambdaH'), info.lambdaH = dtInfo.lambdaH; end
    if isfield(dtInfo,'ifLogRangeEstimate'), info.ifLogRangeEstimate = dtInfo.ifLogRangeEstimate; end
    if isfield(dtInfo,'numRetries'), info.numRetries = dtInfo.numRetries; end
end
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end

function v = local_get_field(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
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

function info = local_small_eps_monitor(W, u, uOld, dt, op, par, ampInfo, gaugeInfo, postGaugeInfo)
info = struct();
if nargin < 7, ampInfo = struct(); end %#ok<NASGU>
if nargin < 8, gaugeInfo = struct(); end %#ok<NASGU>
if nargin < 9, postGaugeInfo = struct(); end %#ok<NASGU>

duNeighbor = diff([u(:); u(1)]);
info.neighborJumpOverEps = max(abs(duNeighbor)) / max(par.eps, realmin);
[dm, dp] = src.one_sided_derivatives(u, op.dtheta, 'first-order');
charSpeed = 2*max(abs([dm(:); dp(:)]));
info.phaseCFL = dt / op.dtheta * charSpeed;
Wpos = W(isfinite(W) & W > 0);
if isempty(Wpos)
    info.logWRange = Inf;
    info.W_max = max(W(:));
else
    info.logWRange = max(log(Wpos(:))) - min(log(Wpos(:)));
    info.W_max = max(Wpos(:));
end
info.uGaugeDrift = max(u(:)) - min(u(:));
% 峰宽/网格比：小于 1--2 时，密度峰形已经无法由当前 theta 网格解析，
% 但 WKB phase 的 argmax 仍可能可靠。
[~, kstar] = max(u(:));
ddu2 = op.Ltheta * u(:);
curv = -real(ddu2(kstar));
if curv > 0
    info.peakWidth = sqrt(par.eps / curv);
    info.peakWidthOverGrid = info.peakWidth / op.dtheta;
else
    info.peakWidth = Inf;
    info.peakWidthOverGrid = Inf;
end
if isstruct(ampInfo) && isfield(ampInfo,'qLogAbsMax')
    info.qLogAbsMax = ampInfo.qLogAbsMax;
elseif isstruct(ampInfo) && isfield(ampInfo,'logQMax') && isfield(ampInfo,'logQMin')
    info.qLogAbsMax = max(abs([ampInfo.logQMax, ampInfo.logQMin]));
    info.logQMax = ampInfo.logQMax;
    info.logQMin = ampInfo.logQMin;
end
if isstruct(gaugeInfo) && isfield(gaugeInfo,'phaseGaugeShift')
    info.phaseGaugeShift = gaugeInfo.phaseGaugeShift;
end
if isstruct(postGaugeInfo) && isfield(postGaugeInfo,'postGaugeShift')
    info.postGaugeShift = postGaugeInfo.postGaugeShift;
end
end

function tf = local_get_bool(s, name, defaultValue)
tf = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end
