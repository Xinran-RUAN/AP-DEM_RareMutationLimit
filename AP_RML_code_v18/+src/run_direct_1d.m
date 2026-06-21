function result = run_direct_1d(par)
%RUN_DIRECT_1D 运行原始密度 n 方程的一维直接求解器。
%   该函数主要用于与 WKB 重构结果对比。因为直接解 n 会在小 eps
%   时出现很窄的 theta 峰，所以在粗 theta 网格上通常更容易欠分辨。

%-------------------------
% 初始化密度、算子、输出目录和进度条
%-------------------------
[n, D, Kx, grid, W0, u0] = model.initial_direct_1d(par);
op = src.build_operators_1d(grid);
par.outdir = utils.ensure_dir(par.outdir);
tag = utils.build_tag(par);
prog = utils.progress_init('DIR-1D', par);
livePlotEverySteps = local_liveplot_every_steps(par);
lp = utils.liveplot_init_1d(par, 'direct');

% 中间结果保存控制器。direct 的 latest checkpoint 使用独立文件，
% 避免覆盖 WKB 的 checkpoint_wkb_latest.mat。
saveCtl = utils.intermediate_save_init(par, [tag '_direct'], 'direct');

history.t = [];
history.residual = [];
history.theta_direct = [];
history.err_direct_to_m = [];
history.sigma_theta = [];
history.mass_total = [];
history.mass_balance_rhs = [];
history.mass_correction_scale = [];
historyEverySteps = local_history_every_steps(par);
lastInfo = struct();
lastPlotState = struct();

t = 0;
step = 0;
residual = Inf;

% 质量平衡修正状态。massState 始终对应最近一个已接受状态。
rhoInitialForMass = op.dtheta * sum(n, 2);
massState = src.mass_balance_stats_1d(rhoInitialForMass, Kx, op);
massInfoCurrent = local_make_initial_mass_info(massState, par);

% 可选保存初始时刻 t=0。若 snapshotTimes 含 0 或 saveInitialSnapshot=true，
% 后处理比较 WKB/direct 时需要 direct 初值也可用。
[initialSaveDue, ~, ~] = utils.intermediate_save_due(saveCtl, t, step);
if initialSaveDue
    try
        rho0 = op.dtheta * sum(n, 2);
        wx0 = utils.trapz_weights(op.nx, op.dx);
        P0 = wx0.' * n;
        [~, idP0] = max(P0);
        theta0 = op.theta(idP0);
        dist0 = utils.periodic_distance(op.theta, theta0);
        diag0 = struct();
        diag0.theta_direct = theta0;
        diag0.err_direct_to_m = utils.periodic_distance(theta0, par.theta_m);
        diag0.sigma_theta = sqrt(sum((dist0.^2).*P0)/max(sum(P0),realmin));
        diag0.rho_min = min(rho0);
        diag0.rho_max = max(rho0);
        diag0.P = P0;
        diag0 = local_attach_mass_info(diag0, massInfoCurrent, massState);
        vars0 = local_make_save_vars(n, rho0, D, Kx, grid, op, par, t, step, residual, history, diag0, tag, W0, u0, massInfoCurrent);
        [saveCtl, didSave0] = utils.intermediate_save_step(saveCtl, par, vars0);
        if didSave0 && par.verbose && isfield(par,'printIntermediateSaveMessage') && par.printIntermediateSaveMessage
            fprintf('[DIR-1D] 已保存初始中间结果: step=0, t=0\n');
        end
    catch ME
        warning('run_direct_1d:InitialSnapshotFailed', '初始 direct 快照保存失败：%s', ME.message);
    end
end

%-------------------------
% 主时间循环：半隐式更新原始 n 方程
%-------------------------
while t < par.T && step < par.maxSteps
    nOld = n;
    nTrial = src.update_direct_density_fd_1d(n, D, Kx, par.eps, par.dt, op, par);

    % 每一步可选进行质量平衡投影。direct 情形下直接缩放 n，rho 同步缩放。
    rhoTrial = op.dtheta * sum(nTrial, 2);
    [massScale, massInfoCurrent] = src.mass_balance_correction_1d( ...
        rhoTrial, Kx, op, par.eps, par.dt, ...
        massState.mass, massState.reactionIntegral, par);
    if massInfoCurrent.applied
        nTrial = massScale * nTrial;
        rhoTrial = massScale * rhoTrial;
    end

    n = nTrial;
    massState = src.mass_balance_stats_1d(rhoTrial, Kx, op);
    massInfoCurrent = local_finalize_mass_info(massInfoCurrent, massState);

    t = t + par.dt;
    step = step + 1;
    residual = max(abs(n(:)-nOld(:))) / (par.dt*(1 + max(abs(n(:)))));

    % 记录 trait marginal 的峰值位置和宽度，并按输出时刻保存完整数据。
    historyDue = (step == 1 || mod(step, historyEverySteps) == 0 || residual <= par.tol);
    plotDue = lp.enabled && (step == 1 || mod(step, livePlotEverySteps) == 0 || residual <= par.tol);
    [saveDue, ~, ~] = utils.intermediate_save_due(saveCtl, t, step);

    if historyDue || plotDue || saveDue
        wx = utils.trapz_weights(op.nx, op.dx);
        P = wx.' * n;
        [~, idP] = max(P);
        theta_direct = op.theta(idP);
        dist = utils.periodic_distance(op.theta, theta_direct);
        sigma_theta = sqrt(sum((dist.^2).*P)/max(sum(P),realmin));
        rhoNow = op.dtheta * sum(n, 2);

        diagNow = struct();
        diagNow.theta_direct = theta_direct;
        diagNow.err_direct_to_m = utils.periodic_distance(theta_direct, par.theta_m);
        diagNow.sigma_theta = sigma_theta;
        diagNow.rho_min = min(rhoNow);
        diagNow.rho_max = max(rhoNow);
        diagNow.P = P;
        diagNow = local_attach_mass_info(diagNow, massInfoCurrent, massState);

        if historyDue
            history.t(end+1,1) = t;
            history.residual(end+1,1) = residual;
            history.theta_direct(end+1,1) = theta_direct;
            history.err_direct_to_m(end+1,1) = diagNow.err_direct_to_m;
            history.sigma_theta(end+1,1) = sigma_theta;
            history.mass_total(end+1,1) = diagNow.mass_total;
            history.mass_correction_scale(end+1,1) = diagNow.mass_correction_scale;
            history.mass_balance_rhs(end+1,1) = diagNow.mass_balance_rhs;
        end

        lastInfo.theta_direct = theta_direct;
        lastInfo.err_direct_to_m = diagNow.err_direct_to_m;
        lastInfo.sigma_theta = sigma_theta;
        lastPlotState = local_make_plot_state(n, rhoNow, Kx, op, par, t, step, residual, P, theta_direct, sigma_theta);

        if saveDue
            vars = local_make_save_vars(n, rhoNow, D, Kx, grid, op, par, t, step, residual, history, diagNow, tag, W0, u0, massInfoCurrent);
            [saveCtl, didSaveIntermediate] = utils.intermediate_save_step(saveCtl, par, vars);
            if didSaveIntermediate && par.verbose && isfield(par,'printIntermediateSaveMessage') && par.printIntermediateSaveMessage
                fprintf('[DIR-1D] 已保存中间结果: step=%d, t=%.6g\n', step, t);
            end
        end
    end

    prog = utils.progress_update(prog, step, t, residual, lastInfo);

    if plotDue && ~isempty(fieldnames(lastPlotState))
        lastPlotState.saveFrame = isfield(par, 'livePlotSave') && par.livePlotSave;
        lp = utils.liveplot_update_1d(lp, 'direct', lastPlotState);
        if isfield(par, 'livePlotPause') && par.livePlotPause > 0
            pause(par.livePlotPause);
        end
    end

    if par.stopByResidual && residual <= par.tol
        break;
    end
end

%-------------------------
% 最终诊断和保存
%-------------------------
rho = op.dtheta * sum(n, 2);
wx = utils.trapz_weights(op.nx, op.dx);
P = wx.' * n;
[~, idP] = max(P);
diagInfo.theta_direct = op.theta(idP);
diagInfo.err_direct_to_m = utils.periodic_distance(diagInfo.theta_direct, par.theta_m);
distFinal = utils.periodic_distance(op.theta, diagInfo.theta_direct);
diagInfo.sigma_theta = sqrt(sum((distFinal.^2).*P)/max(sum(P),realmin));
diagInfo.rho_min = min(rho);
diagInfo.rho_max = max(rho);
diagInfo.P = P;
massState = src.mass_balance_stats_1d(rho, Kx, op);
massInfoCurrent = local_finalize_mass_info(massInfoCurrent, massState);
diagInfo = local_attach_mass_info(diagInfo, massInfoCurrent, massState);
result = struct('n',n,'rho',rho,'D',D,'Kx',Kx,'grid',grid,'op',op,'par',par, ...
                't',t,'step',step,'residual',residual,'history',history, ...
                'diagnostics',diagInfo,'tag',tag,'W0',W0,'u0',u0, ...
                'massInfo',massInfoCurrent);

% 无论是否正常到达终止时刻，最后都补写一次 direct latest checkpoint。
try
    finalVars = local_make_save_vars(n, rho, D, Kx, grid, op, par, t, step, residual, history, diagInfo, tag, W0, u0, massInfoCurrent);
    [saveCtl, ~] = utils.intermediate_save_force_latest(saveCtl, par, finalVars, 'final-or-failed-checkpoint'); %#ok<NASGU>
catch ME
    warning('run_direct_1d:FinalCheckpointFailed', '最终 direct checkpoint 保存失败：%s', ME.message);
end

finishInfo.theta_direct = diagInfo.theta_direct;
finishInfo.err_direct_to_m = diagInfo.err_direct_to_m;
finishInfo.rho_max = diagInfo.rho_max;
utils.progress_finish(prog, step, t, residual, finishInfo);
if local_get_bool(par, 'saveFinalResult', true)
    outfile = utils.solver_file(par.outdir, 'direct', 'final');
    save(outfile, '-struct', 'result');
    if par.verbose
        fprintf('[DIR-1D] 结果已保存: %s\n', outfile);
    end
else
    if par.verbose
        fprintf('[DIR-1D] saveFinalResult=false，仅返回结果并保留已写出的 snapshot。\n');
    end
end
utils.liveplot_finish_1d(lp);
end



function vars = local_make_save_vars(n, rho, D, Kx, grid, op, par, t, step, residual, history, diagInfo, tag, W0, u0, massInfo)
% 组织 direct solver 的完整中间输出变量。
vars = struct();
vars.n = n;
vars.rho = rho;
vars.D = D;
vars.Kx = Kx;
vars.grid = grid;
vars.op = op;
vars.par = par;
vars.t = t;
vars.step = step;
vars.residual = residual;
vars.history = history;
vars.diagnostics = diagInfo;
vars.tag = tag;
vars.W0 = W0;
vars.u0 = u0;
if nargin >= 16 && isstruct(massInfo)
    vars.massInfo = massInfo;
end
vars.createdBy = 'run_direct_1d intermediate output';
end

function info = local_make_initial_mass_info(massState, par)
info = struct();
info.enabled = local_get_bool(par, 'massCorrectionDuringSolve', false);
if isfield(par, 'massCorrection') && ~isempty(par.massCorrection)
    mode = lower(char(par.massCorrection));
    if any(strcmp(mode, {'none','off','false','0','no'})), info.enabled = false; end
    if any(strcmp(mode, {'on','true','1','yes','balance','balance-be','balance-trapezoid'})), info.enabled = true; end
end
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
    diag.mass_correction_scale = local_get_field_num(massInfo, 'scale', 1);
    diag.mass_trial = local_get_field_num(massInfo, 'trialMass', NaN);
    diag.mass_corrected = local_get_field_num(massInfo, 'correctedMass', massState.mass);
else
    diag.mass_correction_enabled = false;
    diag.mass_correction_applied = false;
    diag.mass_correction_scale = 1;
    diag.mass_trial = NaN;
    diag.mass_corrected = massState.mass;
end
end

function v = local_get_field_num(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
end

function livePlotEverySteps = local_liveplot_every_steps(par)
% 根据 par.livePlotEverySteps 或 par.livePlotEveryTime 计算在线图刷新间隔。
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

function state = local_make_plot_state(n, rho, Kx, op, par, t, step, residual, P, theta_direct, sigma_theta)
% 组织 direct solver 在线图像监测所需的数据。
state = struct();
state.n = n;
state.rho = rho;
state.Kx = Kx;
state.x = op.x;
state.theta = op.theta;
state.P = P;
state.theta_peak = theta_direct;
state.theta_m = par.theta_m;
state.sigma_theta = sigma_theta;
state.gammaR = NaN;
state.rho_max = max(rho);
state.residual = residual;
state.t = t;
state.T = par.T;
state.step = step;
state.saveFrame = false;
end

function historyEverySteps = local_history_every_steps(par)
% 根据 par.historyEverySteps 或 par.historyEveryTime 计算诊断记录间隔。
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


function tf = local_get_bool(s, name, defaultValue)
    tf = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        tf = logical(s.(name));
    end
end
