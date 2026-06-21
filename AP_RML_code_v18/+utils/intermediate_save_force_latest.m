function [saveCtl, didSave] = intermediate_save_force_latest(saveCtl, par, vars, reason)
%INTERMEDIATE_SAVE_FORCE_LATEST 强制写一次 latest checkpoint 和 history checkpoint。

if nargin < 4 || isempty(reason)
    reason = 'final-or-failed-checkpoint';
end
if ~isfield(vars,'t') || ~isfield(vars,'step')
    error('intermediate_save_force_latest:InvalidVars', 'vars 必须包含 t 和 step 字段。');
end

didSave = false;
if isfield(saveCtl,'checkpointEnabled') && saveCtl.checkpointEnabled
    varsLatest = vars;
    varsLatest.saveReason = reason;
    varsLatest.createdAt = datestr(now, 31);
    if ~isfield(varsLatest,'createdBy')
        varsLatest.createdBy = 'utils.intermediate_save_force_latest';
    end
    local_save_struct_mat(saveCtl.latestFile, varsLatest, par);
    saveCtl.lastCheckpointStep = vars.step;
    didSave = true;
end

if isfield(saveCtl,'historyEnabled') && saveCtl.historyEnabled
    histVars = struct();
    histVars.par = local_get_field(vars,'par',par);
    histVars.tag = local_get_field(vars,'tag',saveCtl.tag);
    histVars.t = vars.t;
    histVars.step = vars.step;
    histVars.residual = local_get_field(vars,'residual',NaN);
    histVars.history = local_get_field(vars,'history',struct());
    histVars.diagnostics = local_get_field(vars,'diagnostics',struct());
    histVars.dtInfo = local_get_field(vars,'dtInfo',struct());
    histVars.resInfo = local_get_field(vars,'resInfo',struct());
    histVars.createdAt = datestr(now, 31);
    histVars.createdBy = 'utils.intermediate_save_force_latest history checkpoint';
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
    warning('utils:IntermediateSaveFailed', '保存 latest checkpoint 失败：%s。文件：%s', ME.message, filename);
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
