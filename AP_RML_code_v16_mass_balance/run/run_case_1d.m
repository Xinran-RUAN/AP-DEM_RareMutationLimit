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
if ~isfield(cfg, 'caseName') || isempty(cfg.caseName), cfg.caseName = 'case_1d'; end
if ~isfield(cfg, 'dataRoot') || isempty(cfg.dataRoot), cfg.dataRoot = fullfile(root, 'data'); end
cfg.dataRoot = utils.ensure_dir(cfg.dataRoot);

if ~isfield(cfg, 'solvers') || isempty(cfg.solvers), cfg.solvers = {'wkb'}; end
if ~isfield(cfg, 'profile') || isempty(cfg.profile), cfg.profile = 'baseline'; end

% 基本数值参数
cfg = local_default_field(cfg, 'eps', 1e-2);
cfg = local_default_field(cfg, 'Nx', 40);
cfg = local_default_field(cfg, 'Ntheta', 64);
cfg = local_default_field(cfg, 'dt', 1e-3);
cfg = local_default_field(cfg, 'T', 2);
cfg = local_default_field(cfg, 'tol', 1e-8);
cfg = local_default_field(cfg, 'maxSteps', Inf);

% 初值：默认 flat phase、flat W。
if ~isfield(cfg, 'u0_label') || isempty(cfg.u0_label), cfg.u0_label = 'u_flat'; end
if ~isfield(cfg, 'W0_label') || isempty(cfg.W0_label), cfg.W0_label = 'W_flat'; end
if ~isfield(cfg, 'u0_fun') || isempty(cfg.u0_fun), cfg.u0_fun = @(theta) 0.*theta; end
if ~isfield(cfg, 'W0_fun') || isempty(cfg.W0_fun), cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta)); end

% 算法默认值：保持旧版 WKB/direct 可比行为。
cfg = local_default_field(cfg, 'timeIntegrator', 'frozen');
cfg = local_default_field(cfg, 'phaseHamiltonian', 'weno5');
cfg = local_default_field(cfg, 'amplitudeVariant', 'split');
cfg = local_default_field(cfg, 'rhoReconstruction', 'laplace-hybrid');
cfg = local_default_field(cfg, 'reactionDiscretization', 'implicit');
cfg = local_default_field(cfg, 'residualMode', 'legacy');
cfg = local_default_field(cfg, 'massCorrectionDuringSolve', true);
cfg = local_default_field(cfg, 'massCorrectionFormula', 'trapezoid');
cfg = local_default_field(cfg, 'massCorrectionMinScale', 0);
cfg = local_default_field(cfg, 'massCorrectionMaxScale', Inf);
cfg = local_default_field(cfg, 'massCorrectionVerbose', false);
cfg = local_default_field(cfg, 'adaptiveTimeStep', false);
cfg = local_default_field(cfg, 'stopByResidual', false);

% 输出和保存默认值。
cfg = local_default_field(cfg, 'storeSnapshots', true);
cfg = local_default_field(cfg, 'saveIntermediate', true);
cfg = local_default_field(cfg, 'snapshotTimes', []);
cfg = local_default_field(cfg, 'snapshotEveryTime', Inf);
cfg = local_default_field(cfg, 'snapshotEverySteps', []);
cfg = local_default_field(cfg, 'saveInitialSnapshot', true);
cfg = local_default_field(cfg, 'saveInitialAcceptedSnapshot', false);
cfg = local_default_field(cfg, 'saveLatestCheckpoint', true);
cfg = local_default_field(cfg, 'checkpointEveryTime', Inf);
cfg = local_default_field(cfg, 'checkpointEverySteps', []);
cfg = local_default_field(cfg, 'saveHistoryCheckpoint', true);
cfg = local_default_field(cfg, 'historyEveryTime', 0.5);
cfg = local_default_field(cfg, 'historyEverySteps', []);
cfg = local_default_field(cfg, 'saveFinalResult', true);
cfg = local_default_field(cfg, 'saveMatV73', false);

cfg = local_default_field(cfg, 'verbose', true);
cfg = local_default_field(cfg, 'progressEveryPercent', 5);
cfg = local_default_field(cfg, 'progressEverySeconds', 20);
cfg = local_default_field(cfg, 'printIntermediateSaveMessage', false);
cfg = local_default_field(cfg, 'livePlot', false);
cfg = local_default_field(cfg, 'livePlotEveryTime', 0.5);
cfg = local_default_field(cfg, 'livePlotEverySteps', []);
cfg = local_default_field(cfg, 'livePlotSave', false);
cfg = local_default_field(cfg, 'livePlotPause', 0.0);

% 运行控制。
cfg = local_default_field(cfg, 'clearOutput', true);
cfg = local_default_field(cfg, 'continueOnError', false);
cfg = local_default_field(cfg, 'keepResults', false);
cfg = local_default_field(cfg, 'dryRun', false);
end

function cfg = local_default_field(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
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
par = model.default_params_1d(cfg.profile);

% ------------------------- 模型系数 -------------------------
modelFields = {'K_label','D_label','K0','K1','Dmin','D1','theta_m','modelVariant'};
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

% ------------------------- 初值 -------------------------
initialFields = {'u0_label','W0_label','u0_fun','W0_fun', ...
    'initialGaugeMode','prepareInitialMass','initialRhoProfile','initialRhoScale'};
par = local_copy_fields(par, cfg, initialFields);

% ------------------------- 算法选择 -------------------------
algorithmFields = { ...
    'timeIntegrator','phaseHamiltonian','lfAlpha','amplitudeVariant', ...
    'rhoReconstruction','reactionDiscretization','transportImplicit','residualMode', ...
    'massCorrectionDuringSolve','massCorrectionFormula','massCorrectionMinScale','massCorrectionMaxScale','massCorrectionVerbose', ...
    'adaptiveTimeStep','adaptiveStrategy','adaptiveSafety','phaseCflSafety','dtMax','dtMin', ...
    'stopByResidual','smallEpsAutoPatches','enableV13AutoPatches','enableSmallEpsMonitor', ...
    'useSubgridPhasePeak','hamiltonianGauge','ifMode','ifGaugeStrategy','gaugeMode','postGaugeMode', ...
    'thetaRephase','enableStepRetry','maxStepRetries','retryDtFactor','dtMinRetry','rhoLogRetryMax', ...
    'implicitMaxIter','implicitTol','implicitRelax','implicitMinIter','implicitRetryOnFail', ...
    'projectiveMicroSteps','projectiveMicroDtFactor','projectiveMode','projectiveCorrectorSteps'};
par = local_copy_fields(par, cfg, algorithmFields);

if isfield(par, 'dtMax') && isempty(par.dtMax)
    par.dtMax = par.dt;
end

% ------------------------- 保存和输出 -------------------------
saveFields = { ...
    'storeSnapshots','saveIntermediate','snapshotTimes','snapshotEveryTime','snapshotEverySteps', ...
    'saveInitialSnapshot','saveInitialAcceptedSnapshot','saveLatestCheckpoint','checkpointEveryTime','checkpointEverySteps', ...
    'saveHistoryCheckpoint','historyEveryTime','historyEverySteps','saveFinalResult','saveMatV73', ...
    'verbose','progressEveryPercent','progressEverySeconds','progressEverySteps','printIntermediateSaveMessage', ...
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
caseInfo.T = par.T;
caseInfo.profile = par.profile;
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
