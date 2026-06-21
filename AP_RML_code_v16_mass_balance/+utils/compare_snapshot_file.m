function filename = compare_snapshot_file(outdir, index, t)
%COMPARE_SNAPSHOT_FILE Canonical paired comparison snapshot filename.
%  Files are written as:
%      <outdir>/snapshots_compare/snapshot_compare_###_tT.mat
if nargin < 3
    error('utils.compare_snapshot_file:NotEnoughInputs', ...
        'Usage: utils.compare_snapshot_file(outdir, index, t).');
end
index = max(0, round(double(index)));
snapshotDir = utils.ensure_dir(fullfile(char(outdir), 'snapshots_compare'));
filename = fullfile(snapshotDir, sprintf('snapshot_compare_%03d_t%s.mat', index, utils.time_token(t)));
end
