function [saveCtl, didSave] = intermediate_save_step(saveCtl, par, vars)
%INTERMEDIATE_SAVE_STEP 按保存控制器写出中间结果。
%  vars 必须至少包含 t 和 step 字段。其余字段原样保存。

if ~isfield(vars,'t') || ~isfield(vars,'step')
    error('intermediate_save_step:InvalidVars', 'vars 必须包含 t 和 step 字段。');
end

t = vars.t;
step = vars.step;
tolT = 1e-12 + 1e-10*max(1, abs(t));
didSave = false;

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
    varsFull = vars;
    varsFull.saveReason = reason;
    varsFull.createdAt = datestr(now, 31);
    if ~isfield(varsFull,'createdBy')
        varsFull.createdBy = 'utils.intermediate_save_step full snapshot';
    end
    fname = utils.snapshot_filename(saveCtl.outdir, saveCtl.solverName, step, t);
    local_save_struct_mat(fname, varsFull, par);
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
    varsLatest = vars;
    varsLatest.saveReason = 'latest-checkpoint';
    varsLatest.createdAt = datestr(now, 31);
    if ~isfield(varsLatest,'createdBy')
        varsLatest.createdBy = 'utils.intermediate_save_step latest checkpoint';
    end
    local_save_struct_mat(saveCtl.latestFile, varsLatest, par);
    saveCtl.lastCheckpointStep = step;
    didSave = true;
end

if saveCtl.historyEnabled && (saveFull || saveLatest)
    histVars = struct();
    histVars.par = local_get_field(vars,'par',par);
    histVars.tag = local_get_field(vars,'tag',saveCtl.tag);
    histVars.t = t;
    histVars.step = step;
    histVars.residual = local_get_field(vars,'residual',NaN);
    histVars.history = local_get_field(vars,'history',struct());
    histVars.diagnostics = local_get_field(vars,'diagnostics',struct());
    histVars.dtInfo = local_get_field(vars,'dtInfo',struct());
    histVars.resInfo = local_get_field(vars,'resInfo',struct());
    histVars.createdAt = datestr(now, 31);
    histVars.createdBy = 'utils.intermediate_save_step history checkpoint';
    local_save_struct_mat(saveCtl.historyFile, histVars, par);
end
end

function local_save_struct_mat(filename, vars, par)
try
    if local_get_bool(par,'saveMatV73',false)
        save(filename, '-struct', 'vars', '-v7.3');
    else
        save(filename, '-struct', 'vars');
    end
catch ME
    warning('utils:IntermediateSaveFailed', '保存中间结果失败：%s。文件：%s', ME.message, filename);
end
end

function val = local_get_bool(s, field, defaultVal)
if isfield(s, field) && ~isempty(s.(field))
    val = logical(s.(field));
else
    val = defaultVal;
end
end

function val = local_get_field(s, field, defaultVal)
if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
    val = s.(field);
else
    val = defaultVal;
end
end
