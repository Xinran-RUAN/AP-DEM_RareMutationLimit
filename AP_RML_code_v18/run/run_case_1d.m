function runInfo = run_case_1d(cfg)
%RUN_CASE_1D 一维算例统一运行入口。
%
% 用法：
%   cfg = struct();
%   cfg.caseName = 'my_case';
%   cfg.solvers  = {'wkb','direct'};
%   cfg.eps = 1e-2;
%   cfg.Nx = 40;
%   cfg.Ntheta_wkb = 16;
%   cfg.Ntheta_direct = [16, 32];   % 若需要批量跑，把唯一批量变量写成列表
%   cfg.dt = 1e-3;
%   cfg.T  = 2;
%   runInfo = run_case_1d(cfg);
%
% 数据目录固定写成：
%   data/算例名/wkb_eps..._dt..._Nx..._Ntheta..._t...
%   data/算例名/direct_eps..._dt..._Nx..._Ntheta..._t...
%
% 约定：一个 cfg 中最多只有一个真正批量变量；solver 列表不算批量变量。

if nargin < 1 || isempty(cfg)
    cfg = struct();
end
if ~isstruct(cfg)
    error('run_case_1d:InvalidInput', 'cfg 必须是 struct。');
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

cfg = local_fill_defaults(cfg, root);
solvers = local_normalize_solver_list(cfg.solvers);
batchName = local_detect_batch_variable(cfg, solvers);

caseNameForPath = local_safe_token(cfg.caseName, 'case_1d');
caseDir = utils.ensure_dir(fullfile(cfg.dataRoot, caseNameForPath));

runInfo = struct();
runInfo.caseName = cfg.caseName;
runInfo.caseDir = caseDir;
runInfo.dataRoot = cfg.dataRoot;
runInfo.solvers = solvers;
runInfo.batchVariable = batchName;
runInfo.startedAt = datestr(now, 31);
runInfo.finishedAt = '';
runInfo.tasks = struct('solver', {}, 'outdir', {}, 'status', {}, ...
    'eps', {}, 'dt', {}, 'Nx', {}, 'Ntheta', {}, 'T', {}, ...
    'runtimeSeconds', {}, 't', {}, 'step', {}, 'residual', {}, ...
    'resultFile', {}, 'checkpointFile', {}, 'historyFile', {}, ...
    'numSnapshots', {}, 'errorMessage', {});

try
    save(fullfile(caseDir, 'case_config.mat'), 'cfg');
catch ME
    warning('run_case_1d:SaveConfigFailed', 'case_config.mat 保存失败：%s', ME.message);
end

fprintf('\n========== RUN_CASE_1D ==========%s', newline);
fprintf('算例名     : %s%s', cfg.caseName, newline);
fprintf('数据目录   : %s%s', caseDir, newline);
fprintf('求解器     : %s%s', strjoin(solvers, ', '), newline);
if isempty(batchName)
    fprintf('批量变量   : 无%s', newline);
else
    fprintf('批量变量   : %s%s', batchName, newline);
end
fprintf('清空旧目录 : %d%s', logical(cfg.clearOutput), newline);
fprintf('=================================%s', newline);

if logical(cfg.dryRun)
    fprintf('dryRun=true：只打印任务，不真正计算。%s', newline);
end

for is = 1:numel(solvers)
    solverName = solvers{is};
    sweepValues = local_values_for_solver(cfg, solverName, batchName);

    for iv = 1:numel(sweepValues)
        par = local_make_par(cfg, solverName, batchName, sweepValues(iv));
        par.caseName = cfg.caseName;
        par.solverName = solverName;
        par.outdir = local_output_dir(cfg, solverName, par);

        task = local_make_empty_task(solverName, par);
        task.status = 'pending';

        fprintf('\n--- [%s] eps=%.4g, dt=%.4g, Nx=%d, Ntheta=%d, T=%.4g ---%s', ...
            upper(solverName), par.eps, par.dt, par.Nx, par.Ntheta, par.T, newline);
        fprintf('输出目录: %s%s', par.outdir, newline);

        try
            if logical(cfg.clearOutput) && exist(par.outdir, 'dir') == 7
                rmdir(par.outdir, 's');
            end
            utils.ensure_dir(par.outdir);

            caseInfo = local_case_info(cfg, solverName, par);
            save(fullfile(par.outdir, 'case_info.mat'), 'caseInfo', 'par');

            if logical(cfg.dryRun)
                result = [];
                runtimeSeconds = 0;
                task.status = 'dry-run';
            else
                ticRun = tic;
                switch solverName
                    case 'wkb'
                        result = src.run_wkb_1d(par);
                    case 'direct'
                        result = src.run_direct_1d(par);
                    otherwise
                        error('run_case_1d:UnknownSolver', '未知求解器：%s', solverName);
                end
                runtimeSeconds = toc(ticRun);
                task.status = 'ok';
                task = local_fill_task_from_result(task, result);
            end

            task.runtimeSeconds = runtimeSeconds;
            task.resultFile = local_existing_or_empty(utils.solver_file(par.outdir, solverName, 'final'));
            task.checkpointFile = local_existing_or_empty(utils.solver_file(par.outdir, solverName, 'checkpoint'));
            task.historyFile = local_existing_or_empty(utils.solver_file(par.outdir, solverName, 'history'));
            task.numSnapshots = local_count_files(par.outdir, sprintf('snapshot_%s_*.mat', solverName));

            if logical(cfg.keepResults)
                runInfo.results(is, iv).solver = solverName; %#ok<AGROW>
                runInfo.results(is, iv).data = result; %#ok<AGROW>
            end

            fprintf('[完成] %s: runtime=%.2f s, snapshots=%d%s', ...
                solverName, task.runtimeSeconds, task.numSnapshots, newline);

        catch ME
            task.status = 'failed';
            task.errorMessage = ME.message;
            fprintf(2, '[失败] %s: %s%s', solverName, ME.message, newline);

            runInfo.tasks(end+1) = task; %#ok<AGROW>
            local_save_run_summary(caseDir, runInfo, cfg);

            if ~logical(cfg.continueOnError)
                rethrow(ME);
            end
            continue;
        end

        runInfo.tasks(end+1) = task; %#ok<AGROW>
        local_save_run_summary(caseDir, runInfo, cfg);
    end
end

runInfo.finishedAt = datestr(now, 31);
local_save_run_summary(caseDir, runInfo, cfg);

fprintf('\nRUN_CASE_1D 完成。summary: %s%s', fullfile(caseDir, 'run_summary.mat'), newline);
end

% =====================================================================
% 参数整理
% =====================================================================

function cfg = local_fill_defaults(cfg, root)
% 本工程不再使用 cfg.profile 或 +model/default_params_1d.m。
% run_wkb_main_1d.m 和 run_direct_main_1d.m 必须显式给出所有用户参数。
% 这里主要做完整性检查，避免参数被隐藏在其它文件中。

if isfield(cfg, 'profile')
    error('run_case_1d:ProfileRemoved', ...
        ['cfg.profile 已删除。请不要使用 baseline/profile 预设；', ...
         '所有 K、D、初值、网格、时间步和算法参数都应写在主程序开头。']);
end

% 为了兼容用户直接从命令行构造 cfg，dataRoot 允许唯一一个路径默认。
% 两个主程序中已经显式写了 cfg.dataRoot。
if ~isfield(cfg, 'dataRoot') || isempty(cfg.dataRoot)
    cfg.dataRoot = fullfile(root, 'data');
end
cfg.dataRoot = utils.ensure_dir(cfg.dataRoot);

basicFields = { ...
    'caseName','solvers','eps','Nx','dt','T','tol','maxSteps', ...
    'K_label','D_label','K0','K1','Dmin','D1','theta_m','K_fun','D_fun', ...
    'u0_label','W0_label','u0_fun','W0_fun','initialGaugeMode','prepareInitialMass', ...
    'massCorrectionDuringSolve','massCorrectionFormula','massCorrectionMinScale','massCorrectionMaxScale','massCorrectionVerbose', ...
    'adaptiveTimeStep','stopByResidual','progressPrintWhenTolReached', ...
    'storeSnapshots','saveIntermediate','snapshotTimes','snapshotEveryTime','snapshotEverySteps', ...
    'saveInitialSnapshot','saveInitialAcceptedSnapshot','saveLatestCheckpoint','checkpointEveryTime','checkpointEverySteps', ...
    'saveHistoryCheckpoint','historyEveryTime','historyEverySteps','saveFinalResult','saveMatV73', ...
    'clearOutput','continueOnError','keepResults','dryRun', ...
    'verbose','progressEveryPercent','progressEverySeconds','progressEverySteps','progressPrintWhenTolReached','printIntermediateSaveMessage', ...
    'livePlot','livePlotEveryTime','livePlotEverySteps','livePlotSave','livePlotPause'};
local_require_fields(cfg, basicFields, '通用参数');

solvers = local_normalize_solver_list(cfg.solvers);
if any(strcmp(solvers, 'wkb'))
    wkbFields = { ...
        'Ntheta_wkb','useDualThetaGrid','Ntheta_u','Ntheta_W', ...
        'timeIntegrator','phaseHamiltonian','lfAlpha','amplitudeVariant','transportImplicit', ...
        'rhoReconstruction','reactionDiscretization','residualMode', ...
        'clipNegativeW','WFloor','monitorMatrixCondition', ...
        'HSolverMode','HDtol','HInterpND','HLocalRadius','eigenSolver','eigenDenseThreshold','eigenTol','eigenMaxit', ...
        'dualHSolverModeW','hamiltonianGauge','fullEigenRelaxHMode', ...
        'dualWtoUInterp','dualUtoWInterp','dualHtoWInterp', ...
        'adaptiveStrategy','adaptiveSafety','phaseCflSafety','dtMax','dtMin'};
    local_require_fields(cfg, wkbFields, 'WKB 参数');
end
if any(strcmp(solvers, 'direct'))
    directFields = {'Ntheta_direct','reactionDiscretization','residualMode'};
    local_require_fields(cfg, directFields, 'direct 参数');
end

% 不再允许 2D/旧版预设相关字段混入主程序。
forbidden = {'modelVariant','enableV13AutoPatches','smallEpsAutoPatches','patchPreset'};
for i = 1:numel(forbidden)
    if isfield(cfg, forbidden{i})
        error('run_case_1d:ForbiddenPresetField', ...
            '字段 cfg.%s 属于旧版预设/补丁路径，已删除。请在主程序中显式设置需要的参数。', forbidden{i});
    end
end
end

function local_require_fields(cfg, fields, groupName)
missing = {};
for i = 1:numel(fields)
    f = fields{i};
    if ~isfield(cfg, f) || isempty(cfg.(f))
        % [] 作为显式参数时，有少数字段是合法的。
        allowEmpty = any(strcmp(f, {'Ntheta_W','snapshotEverySteps','checkpointEverySteps','historyEverySteps','progressEverySteps','livePlotEverySteps'}));
        if ~(isfield(cfg, f) && allowEmpty)
            missing{end+1} = f; %#ok<AGROW>
        end
    end
end
if ~isempty(missing)
    error('run_case_1d:MissingExplicitParameters', ...
        '%s缺少显式参数：%s。请在 run_wkb_main_1d.m 或 run_direct_main_1d.m 开头写明。', ...
        groupName, strjoin(missing, ', '));
end
end

function solvers = local_normalize_solver_list(solvers)
if ischar(solvers) || local_isstring(solvers)
    solvers = cellstr(solvers);
end
if ~iscell(solvers) || isempty(solvers)
    error('run_case_1d:InvalidSolvers', 'cfg.solvers 必须是字符串或 cell。');
end
for i = 1:numel(solvers)
    solvers{i} = utils.normalize_solver_name(solvers{i});
    if ~any(strcmp(solvers{i}, {'wkb','direct'}))
        error('run_case_1d:InvalidSolver', '目前 run_case_1d 只支持 wkb/direct，收到：%s', solvers{i});
    end
end
solvers = unique(solvers, 'stable');
end

function batchName = local_detect_batch_variable(cfg, solvers)
% 一个 cfg 中最多允许一个真正的数值批量变量。
% solvers={'wkb','direct'} 不算批量变量。
candidates = {};
commonFields = {'eps','Nx','dt','T'};
for i = 1:numel(commonFields)
    f = commonFields{i};
    if isfield(cfg, f) && isnumeric(cfg.(f)) && numel(cfg.(f)) > 1
        candidates{end+1} = f; %#ok<AGROW>
    end
end

hasNthetaBatch = false;
if isfield(cfg, 'Ntheta') && isnumeric(cfg.Ntheta) && numel(cfg.Ntheta) > 1
    hasNthetaBatch = true;
end
for i = 1:numel(solvers)
    f = ['Ntheta_' solvers{i}];
    if isfield(cfg, f) && isnumeric(cfg.(f)) && numel(cfg.(f)) > 1
        hasNthetaBatch = true;
    end
end
if hasNthetaBatch
    candidates{end+1} = 'Ntheta'; %#ok<AGROW>
end

candidates = unique(candidates, 'stable');
if numel(candidates) > 1
    error('run_case_1d:TooManyBatchVariables', ...
        '一个主程序最多只批量跑一个变量。当前批量变量为：%s。', strjoin(candidates, ', '));
end
if isempty(candidates)
    batchName = '';
else
    batchName = candidates{1};
end
end

function values = local_values_for_solver(cfg, solverName, batchName)
if isempty(batchName)
    values = 1;
    return;
end
if strcmp(batchName, 'Ntheta')
    values = local_ntheta_values(cfg, solverName);
else
    values = cfg.(batchName)(:).';
end
end

function values = local_ntheta_values(cfg, solverName)
fieldName = ['Ntheta_' solverName];
if isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
    values = cfg.(fieldName)(:).';
else
    values = cfg.Ntheta(:).';
end
end

function par = local_make_par(cfg, solverName, batchName, batchValue)
par = struct();

% ------------------------- 模型系数 -------------------------
modelFields = {'K_label','D_label','K0','K1','Dmin','D1','theta_m'};
par = local_copy_fields(par, cfg, modelFields);

if isfield(cfg, 'K_fun') && ~isempty(cfg.K_fun)
    par.K_fun = cfg.K_fun;
elseif isfield(par, 'K0') && isfield(par, 'K1')
    K0 = par.K0;
    K1 = par.K1;
    par.K_fun = @(x) K0 - K1*cos(2*pi*x);
end

if isfield(cfg, 'D_fun') && ~isempty(cfg.D_fun)
    par.D_fun = cfg.D_fun;
elseif isfield(par, 'Dmin') && isfield(par, 'D1') && isfield(par, 'theta_m')
    Dmin = par.Dmin;
    D1 = par.D1;
    theta_m = par.theta_m;
    par.D_fun = @(theta) Dmin + D1*(1 - cos(2*pi*(theta - theta_m)));
end

% ------------------------- 基本参数 -------------------------
par.eps = local_scalar_or_batch(cfg, 'eps', batchName, batchValue);
par.Nx  = round(local_scalar_or_batch(cfg, 'Nx', batchName, batchValue));
par.dt  = local_scalar_or_batch(cfg, 'dt', batchName, batchValue);
par.T   = local_scalar_or_batch(cfg, 'T', batchName, batchValue);
par.tol = local_first_scalar(cfg.tol);
par.maxSteps = local_first_scalar(cfg.maxSteps);

if strcmp(batchName, 'Ntheta')
    par.Ntheta = round(batchValue);
else
    nthetaValues = local_ntheta_values(cfg, solverName);
    if numel(nthetaValues) ~= 1
        error('run_case_1d:NthetaListWithoutBatch', ...
            'Ntheta 是列表时会作为唯一批量变量；请确认没有其它批量变量。');
    end
    par.Ntheta = round(nthetaValues(1));
end

% 双网格约定：par.Ntheta 是 W 的粗网格点数，用于目录名和 W 更新；
% par.Ntheta_u 是 phase 的细网格点数。只对 WKB 生效。
if strcmp(solverName, 'wkb') && isfield(cfg, 'useDualThetaGrid') && logical(cfg.useDualThetaGrid)
    par.useDualThetaGrid = true;
    if isfield(cfg, 'Ntheta_W') && ~isempty(cfg.Ntheta_W)
        par.Ntheta_W = round(local_first_scalar(cfg.Ntheta_W));
        par.Ntheta = par.Ntheta_W;
    else
        par.Ntheta_W = par.Ntheta;
    end
    if isfield(cfg, 'Ntheta_u') && ~isempty(cfg.Ntheta_u)
        par.Ntheta_u = round(local_first_scalar(cfg.Ntheta_u));
    else
        par.Ntheta_u = max(par.Ntheta_W, 4*par.Ntheta_W);
    end
else
    par.useDualThetaGrid = false;
end

% ------------------------- 初值 -------------------------
initialFields = {'u0_label','W0_label','u0_fun','W0_fun', ...
    'initialGaugeMode','prepareInitialMass','initialRhoProfile','initialRhoScale'};
par = local_copy_fields(par, cfg, initialFields);

% ------------------------- 算法选择 -------------------------
algorithmFields = { ...
    'timeIntegrator','phaseHamiltonian','lfAlpha','amplitudeVariant', ...
    'rhoReconstruction','reactionDiscretization','transportImplicit','residualMode', ...
    'clipNegativeW','WFloor','monitorMatrixCondition', ...
    'massCorrectionDuringSolve','massCorrectionFormula','massCorrectionMinScale','massCorrectionMaxScale','massCorrectionVerbose', ...
    'HSolverMode','HDtol','HInterpND','HLocalRadius','eigenSolver','eigenDenseThreshold','eigenTol','eigenMaxit', ...
    'useDualThetaGrid','Ntheta_u','Ntheta_W','dualWtoUInterp','dualUtoWInterp','dualHtoWInterp','dualHSolverModeW','fullEigenRelaxHMode', ...
    'adaptiveTimeStep','adaptiveStrategy','adaptiveSafety','phaseCflSafety','dtMax','dtMin', ...
    'stopByResidual','hamiltonianGauge'};
par = local_copy_fields(par, cfg, algorithmFields);

if isfield(par, 'dtMax') && isempty(par.dtMax)
    par.dtMax = par.dt;
end

% ------------------------- 保存和输出 -------------------------
saveFields = { ...
    'storeSnapshots','saveIntermediate','snapshotTimes','snapshotEveryTime','snapshotEverySteps', ...
    'saveInitialSnapshot','saveInitialAcceptedSnapshot','saveLatestCheckpoint','checkpointEveryTime','checkpointEverySteps', ...
    'saveHistoryCheckpoint','historyEveryTime','historyEverySteps','saveFinalResult','saveMatV73', ...
    'verbose','progressEveryPercent','progressEverySeconds','progressEverySteps','progressPrintWhenTolReached','printIntermediateSaveMessage', ...
    'livePlot','livePlotEveryTime','livePlotEverySteps','livePlotSave','livePlotPause'};
par = local_copy_fields(par, cfg, saveFields);

% snapshotTimes 只保留 [0,T] 内的时刻。
if isfield(par, 'snapshotTimes') && ~isempty(par.snapshotTimes)
    par.snapshotTimes = unique(par.snapshotTimes(:).');
    par.snapshotTimes = par.snapshotTimes(par.snapshotTimes >= -1e-12 & par.snapshotTimes <= par.T + 1e-12);
end

% direct 不需要 time-AP，但保留字段便于 case_info 中追溯。
par.requestedSolver = solverName;
end

function par = local_copy_fields(par, cfg, fields)
for i = 1:numel(fields)
    f = fields{i};
    if isfield(cfg, f) && ~isempty(cfg.(f))
        par.(f) = cfg.(f);
    end
end
end

function v = local_scalar_or_batch(cfg, name, batchName, batchValue)
if strcmp(batchName, name)
    v = batchValue;
    return;
end
v = local_first_scalar(cfg.(name));
end

function v = local_first_scalar(x)
if isempty(x)
    v = [];
    return;
end
if isnumeric(x)
    v = x(1);
else
    v = x;
end
end

function outdir = local_output_dir(cfg, solverName, par)
caseNameForPath = local_safe_token(cfg.caseName, 'case_1d');
caseDir = fullfile(cfg.dataRoot, caseNameForPath);
solverDir = sprintf('%s_eps%s_dt%s_Nx%d_Ntheta%d_t%s', ...
    solverName, local_num_token(par.eps), local_num_token(par.dt), ...
    par.Nx, par.Ntheta, local_num_token(par.T));
if strcmp(solverName, 'wkb') && isfield(par, 'useDualThetaGrid') && logical(par.useDualThetaGrid)
    solverDir = sprintf('%s_Nu%d', solverDir, par.Ntheta_u);
end
if isfield(cfg, 'dirExtraTag') && ~isempty(cfg.dirExtraTag)
    solverDir = [solverDir '_' local_safe_token(cfg.dirExtraTag, 'tag')];
end
outdir = fullfile(caseDir, solverDir);
end

% =====================================================================
% 记录和保存
% =====================================================================

function task = local_make_empty_task(solverName, par)
task = struct();
task.solver = solverName;
task.outdir = par.outdir;
task.status = '';
task.eps = par.eps;
task.dt = par.dt;
task.Nx = par.Nx;
task.Ntheta = par.Ntheta;
task.T = par.T;
task.runtimeSeconds = NaN;
task.t = NaN;
task.step = NaN;
task.residual = NaN;
task.resultFile = '';
task.checkpointFile = '';
task.historyFile = '';
task.numSnapshots = 0;
task.errorMessage = '';
end

function task = local_fill_task_from_result(task, result)
if isstruct(result)
    if isfield(result, 't'), task.t = result.t; end
    if isfield(result, 'step'), task.step = result.step; end
    if isfield(result, 'residual'), task.residual = result.residual; end
end
end

function caseInfo = local_case_info(cfg, solverName, par)
caseInfo = struct();
caseInfo.caseName = cfg.caseName;
caseInfo.solver = solverName;
caseInfo.outdir = par.outdir;
caseInfo.createdAt = datestr(now, 31);
caseInfo.eps = par.eps;
caseInfo.dt = par.dt;
caseInfo.Nx = par.Nx;
caseInfo.Ntheta = par.Ntheta;
caseInfo.useDualThetaGrid = local_get_field(par, 'useDualThetaGrid', false);
caseInfo.Ntheta_u = local_get_field(par, 'Ntheta_u', []);
caseInfo.Ntheta_W = local_get_field(par, 'Ntheta_W', []);
caseInfo.HSolverMode = local_get_text(par, 'HSolverMode', '');
caseInfo.dualHSolverModeW = local_get_text(par, 'dualHSolverModeW', '');
caseInfo.fullEigenRelaxHMode = local_get_text(par, 'fullEigenRelaxHMode', '');
caseInfo.T = par.T;
caseInfo.u0_label = local_get_text(par, 'u0_label', '');
caseInfo.W0_label = local_get_text(par, 'W0_label', '');
caseInfo.K_label = local_get_text(par, 'K_label', '');
caseInfo.D_label = local_get_text(par, 'D_label', '');
caseInfo.timeIntegrator = local_get_text(par, 'timeIntegrator', '');
caseInfo.phaseHamiltonian = local_get_text(par, 'phaseHamiltonian', '');
caseInfo.amplitudeVariant = local_get_text(par, 'amplitudeVariant', '');
caseInfo.rhoReconstruction = local_get_text(par, 'rhoReconstruction', '');
caseInfo.reactionDiscretization = local_get_text(par, 'reactionDiscretization', '');
caseInfo.massCorrectionDuringSolve = local_get_field(par, 'massCorrectionDuringSolve', []);
caseInfo.massCorrectionFormula = local_get_text(par, 'massCorrectionFormula', '');
caseInfo.snapshotTimes = local_get_field(par, 'snapshotTimes', []);
caseInfo.dataLayout = 'data/caseName/solver_eps..._dt..._Nx..._Ntheta..._t...';
end

function local_save_run_summary(caseDir, runInfo, cfg)
try
    save(fullfile(caseDir, 'run_summary.mat'), 'runInfo', 'cfg');
catch ME
    warning('run_case_1d:SaveSummaryFailed', 'run_summary.mat 保存失败：%s', ME.message);
end
end

function path = local_existing_or_empty(path)
if exist(path, 'file') ~= 2
    path = '';
end
end

function n = local_count_files(folder, pattern)
if exist(folder, 'dir') ~= 7
    n = 0;
    return;
end
files = dir(fullfile(folder, pattern));
n = numel(files);
end

function value = local_get_field(s, name, defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end

function txt = local_get_text(s, name, defaultText)
txt = defaultText;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    txt = char(s.(name));
end
end

% =====================================================================
% 文件名 token
% =====================================================================

function s = local_safe_token(v, defaultValue)
if nargin < 2
    defaultValue = 'token';
end
if isempty(v)
    s = defaultValue;
    return;
end
s = local_to_char(v);
s = strtrim(s);
% 只替换常见路径非法字符和空白，中文算例名会保留。
s = regexprep(s, '[\\/:*?"<>|\s]+', '_');
s = regexprep(s, '^_+|_+$', '');
if isempty(s)
    s = defaultValue;
end
end

function tf = local_isstring(v)
tf = false;
if exist('isstring', 'builtin') || exist('isstring', 'file')
    tf = isstring(v);
end
end

function s = local_to_char(v)
if ischar(v)
    s = v;
elseif local_isstring(v)
    s = char(v);
elseif isnumeric(v) && isscalar(v)
    s = num2str(v);
else
    try
        s = char(v);
    catch
        s = 'token';
    end
end
end

function s = local_num_token(x)
% 使用统一的数值命名规则，避免 2.5e-4 被四舍五入成 3e-4。
s = utils.num_token(x);
end
