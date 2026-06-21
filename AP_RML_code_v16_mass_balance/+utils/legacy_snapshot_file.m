function filename = legacy_snapshot_file(outdir, solverName, Nx, Ntheta, t)
%LEGACY_SNAPSHOT_FILE Canonical snapshot filename for legacy 2-D scripts.
%
%  快照文件直接放在 outdir 下，不再创建 snapshots 子文件夹。
%  WKB legacy 2-D snapshots:
%      <outdir>/snapshot_legacy_wkb_K0015_t2.mat
%  Direct legacy snapshots:
%      <outdir>/snapshot_legacy_direct_Nx0005_K0500_t2.mat
%
%  Nx may be [] for mesh-based WKB runs where no 1-D Nx exists.
if nargin < 5
    error('utils.legacy_snapshot_file:NotEnoughInputs', ...
        'Usage: utils.legacy_snapshot_file(outdir, solverName, Nx, Ntheta, t).');
end

solver = utils.normalize_solver_name(solverName);
outdir = utils.ensure_dir(char(outdir));
Ntheta = max(0, round(double(Ntheta)));

if isempty(Nx) || ~isfinite(double(Nx))
    name = sprintf('snapshot_%s_K%04d_t%s.mat', solver, Ntheta, utils.time_token(t));
else
    Nx = max(0, round(double(Nx)));
    name = sprintf('snapshot_%s_Nx%04d_K%04d_t%s.mat', solver, Nx, Ntheta, utils.time_token(t));
end

filename = fullfile(outdir, name);
end
