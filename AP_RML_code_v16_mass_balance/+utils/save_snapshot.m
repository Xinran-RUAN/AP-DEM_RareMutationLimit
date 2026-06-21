function save_snapshot(outdir, tag, t, vars)
%SAVE_SNAPSHOT Save a generic snapshot directly under outdir.
%
%  Generic snapshots are written as:
%      <outdir>/snapshot_<tag>_tT.mat
if nargin < 4 || isempty(vars)
    vars = struct();
end

safeTag = regexprep(char(tag), '[^A-Za-z0-9_\.-]', '_');
outdir = utils.ensure_dir(char(outdir));
filename = fullfile(outdir, sprintf('snapshot_%s_t%s.mat', safeTag, utils.time_token(t)));
save(filename, '-struct', 'vars');
end
