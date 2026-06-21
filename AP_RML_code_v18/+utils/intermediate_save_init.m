function saveCtl = intermediate_save_init(par, tag, solverName)
%INTERMEDIATE_SAVE_INIT 初始化中间结果保存控制器。
%
%  设计目的：
%    1. full snapshot 直接保存到 par.outdir/，统一命名为 snapshot_<solver>_stepXXXXXXX_tT.mat；
%    2. latest checkpoint 定期覆盖保存，统一命名为 checkpoint_<solver>_latest.mat；
%    3. history checkpoint 只保存 history 和轻量诊断，统一命名为 history_<solver>_latest.mat。
%
%  solverName 用于区分 WKB 与 direct 的 latest checkpoint，避免 direct
%  覆盖 WKB 的 checkpoint_wkb_latest.mat。WKB 默认使用 checkpoint_wkb_latest.mat，
%  direct 默认使用 checkpoint_direct_latest.mat。

if nargin < 3 || isempty(solverName)
    solverName = 'wkb';
end
solverName = lower(char(solverName));

saveCtl = struct();
saveCtl.tag = tag;
saveCtl.solverName = solverName;
saveCtl.outdir = utils.ensure_dir(par.outdir);
% 快照文件直接写入 outdir，不再创建 snapshots 子文件夹。
saveCtl.snapshotDir = saveCtl.outdir;

saveCtl.snapshotEnabled = local_get_bool(par,'storeSnapshots',false) || ...
                          local_get_bool(par,'saveIntermediateSnapshots',false) || ...
                          local_get_bool(par,'saveIntermediate',false) || ...
                          local_get_bool(par,'saveInitialSnapshot',false);
if isfield(par,'saveIntermediate') && ~isempty(par.saveIntermediate)
    saveCtl.snapshotEnabled = local_get_bool(par,'saveIntermediate',saveCtl.snapshotEnabled) || local_get_bool(par,'saveInitialSnapshot',false);
end
saveCtl.checkpointEnabled = local_get_bool(par,'saveLatestCheckpoint',true);
saveCtl.historyEnabled = local_get_bool(par,'saveHistoryCheckpoint',true);

if isfield(par,'snapshotTimes') && ~isempty(par.snapshotTimes)
    saveCtl.snapshotTimes = sort(par.snapshotTimes(:).');
else
    saveCtl.snapshotTimes = [];
end
if local_get_bool(par,'saveInitialSnapshot',false)
    saveCtl.snapshotTimes = unique([0, saveCtl.snapshotTimes]);
end
saveCtl.nextSnapshotIndex = 1;

saveCtl.snapshotEveryTime = local_get_num(par,'snapshotEveryTime',Inf);
if ~(isfinite(saveCtl.snapshotEveryTime) && saveCtl.snapshotEveryTime > 0)
    saveCtl.snapshotEveryTime = Inf;
end
saveCtl.snapshotEverySteps = round(local_get_num(par,'snapshotEverySteps',Inf));
if ~(isfinite(saveCtl.snapshotEverySteps) && saveCtl.snapshotEverySteps > 0)
    saveCtl.snapshotEverySteps = Inf;
end
saveCtl.nextPeriodicSnapshotTime = saveCtl.snapshotEveryTime;

saveCtl.checkpointEveryTime = local_get_num(par,'checkpointEveryTime',1.0);
if ~(isfinite(saveCtl.checkpointEveryTime) && saveCtl.checkpointEveryTime > 0)
    saveCtl.checkpointEveryTime = Inf;
end
saveCtl.checkpointEverySteps = round(local_get_num(par,'checkpointEverySteps',Inf));
if ~(isfinite(saveCtl.checkpointEverySteps) && saveCtl.checkpointEverySteps > 0)
    saveCtl.checkpointEverySteps = Inf;
end
saveCtl.nextCheckpointTime = saveCtl.checkpointEveryTime;

switch solverName
    case {'wkb','wkb-1d'}
        defaultCheckpointFile = 'checkpoint_wkb_latest.mat';
        defaultHistoryFile = 'history_wkb_latest.mat';
    otherwise
        defaultCheckpointFile = sprintf('checkpoint_%s_latest.mat', utils.normalize_solver_name(solverName));
        defaultHistoryFile = sprintf('history_%s_latest.mat', utils.normalize_solver_name(solverName));
end
saveCtl.latestFile = fullfile(saveCtl.outdir, local_get_string(par,'checkpointFile',defaultCheckpointFile));
saveCtl.historyFile = fullfile(saveCtl.outdir, local_get_string(par,'historyCheckpointFile',defaultHistoryFile));

% 若是 direct 或其它辅助求解器，默认不要覆盖 WKB 的 checkpoint。
% 只有用户显式设置 par.directCheckpointFile / par.directHistoryCheckpointFile
% 时才使用指定文件名；否则使用 checkpoint_direct_latest.mat 等独立文件。
if ~any(strcmp(solverName, {'wkb','wkb-1d'}))
    solverPrefix = matlab.lang.makeValidName(solverName);
    ckField = [solverPrefix 'CheckpointFile'];
    histField = [solverPrefix 'HistoryCheckpointFile'];
    if isfield(par, ckField) && ~isempty(par.(ckField))
        saveCtl.latestFile = fullfile(saveCtl.outdir, char(par.(ckField)));
    else
        saveCtl.latestFile = fullfile(saveCtl.outdir, defaultCheckpointFile);
    end
    if isfield(par, histField) && ~isempty(par.(histField))
        saveCtl.historyFile = fullfile(saveCtl.outdir, char(par.(histField)));
    else
        saveCtl.historyFile = fullfile(saveCtl.outdir, defaultHistoryFile);
    end
end

saveCtl.lastSnapshotStep = -Inf;
saveCtl.lastCheckpointStep = -Inf;
saveCtl.savedCount = 0;
end

function val = local_get_bool(s, field, defaultVal)
if isfield(s, field) && ~isempty(s.(field))
    val = logical(s.(field));
else
    val = defaultVal;
end
end

function val = local_get_num(s, field, defaultVal)
if isfield(s, field) && ~isempty(s.(field)) && isnumeric(s.(field)) && isscalar(s.(field))
    val = s.(field);
else
    val = defaultVal;
end
end

function val = local_get_string(s, field, defaultVal)
if isfield(s, field) && ~isempty(s.(field))
    val = char(s.(field));
else
    val = defaultVal;
end
end
