function filename = solver_file(outdir, solverName, kind)
%SOLVER_FILE Canonical filenames for solver-level data.
%
%  result files:
%      result_wkb_final.mat
%      result_direct_final.mat
%  latest checkpoints:
%      checkpoint_wkb_latest.mat
%      checkpoint_direct_latest.mat
%  history checkpoints:
%      history_wkb_latest.mat
%      history_direct_latest.mat
if nargin < 3
    error('utils.solver_file:NotEnoughInputs', 'Usage: utils.solver_file(outdir, solverName, kind).');
end
solver = utils.normalize_solver_name(solverName);
kind = lower(strtrim(char(kind)));
switch kind
    case {'final','result','result-final'}
        name = sprintf('result_%s_final.mat', solver);
    case {'checkpoint','latest','latest-checkpoint'}
        name = sprintf('checkpoint_%s_latest.mat', solver);
    case {'history','history-checkpoint'}
        name = sprintf('history_%s_latest.mat', solver);
    otherwise
        error('utils.solver_file:UnknownKind', 'Unknown solver file kind "%s".', kind);
end
filename = fullfile(char(outdir), name);
end
