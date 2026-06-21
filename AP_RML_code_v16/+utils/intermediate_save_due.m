function [needSave, needFull, needLatest] = intermediate_save_due(saveCtl, t, step)
%INTERMEDIATE_SAVE_DUE 判断当前时刻是否需要保存中间结果。
%  只做判断，不更新 saveCtl。真正保存和更新阈值由 intermediate_save_step 完成。

if nargin < 3
    step = 0;
end

tolT = 1e-12 + 1e-10*max(1, abs(t));
needFull = false;
needLatest = false;

if isfield(saveCtl,'snapshotEnabled') && saveCtl.snapshotEnabled
    if step == 1
        needFull = true;
    end
    if isfield(saveCtl,'nextSnapshotIndex') && isfield(saveCtl,'snapshotTimes')
        if saveCtl.nextSnapshotIndex <= numel(saveCtl.snapshotTimes) && ...
                t + tolT >= saveCtl.snapshotTimes(saveCtl.nextSnapshotIndex)
            needFull = true;
        end
    end
    if isfield(saveCtl,'nextPeriodicSnapshotTime') && isfinite(saveCtl.nextPeriodicSnapshotTime) && ...
            t + tolT >= saveCtl.nextPeriodicSnapshotTime
        needFull = true;
    end
    if isfield(saveCtl,'snapshotEverySteps') && isfinite(saveCtl.snapshotEverySteps) && ...
            saveCtl.snapshotEverySteps > 0 && step > 0 && mod(step, saveCtl.snapshotEverySteps) == 0 && ...
            step ~= saveCtl.lastSnapshotStep
        needFull = true;
    end
end

if isfield(saveCtl,'checkpointEnabled') && saveCtl.checkpointEnabled
    if step == 1
        needLatest = true;
    end
    if isfield(saveCtl,'nextCheckpointTime') && isfinite(saveCtl.nextCheckpointTime) && ...
            t + tolT >= saveCtl.nextCheckpointTime
        needLatest = true;
    end
    if isfield(saveCtl,'checkpointEverySteps') && isfinite(saveCtl.checkpointEverySteps) && ...
            saveCtl.checkpointEverySteps > 0 && step > 0 && mod(step, saveCtl.checkpointEverySteps) == 0 && ...
            step ~= saveCtl.lastCheckpointStep
        needLatest = true;
    end
end

needSave = needFull || needLatest;
end
