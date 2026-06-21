function filename = snapshot_filename(outdir, solverName, step, t)
%SNAPSHOT_FILENAME Canonical full snapshot filename.
%
%  快照文件直接放在 outdir 下，不再创建 snapshots 子文件夹。
%  新版文件名只保留 solver 和时间，不再保留 step 编号：
%      <outdir>/snapshot_<solver>_t5.mat
%
%  输入参数 step 仅为兼容旧调用接口而保留，文件名中不使用。
if nargin < 4
    error('utils.snapshot_filename:NotEnoughInputs', ...
        'Usage: utils.snapshot_filename(outdir, solverName, step, t).');
end

% step intentionally unused. Keep it to avoid changing calls elsewhere.
%#ok<NASGU>

solver = utils.normalize_solver_name(solverName);
outdir = utils.ensure_dir(char(outdir));
filename = fullfile(outdir, sprintf('snapshot_%s_t%s.mat', ...
    solver, utils.time_token(t)));
end
