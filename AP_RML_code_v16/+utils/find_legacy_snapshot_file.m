function [fileName, selectedTime, found] = find_legacy_snapshot_file(dataDir, solverName, Nx, Ntheta, targetTime, mustExist)
%FIND_LEGACY_SNAPSHOT_FILE Locate canonical legacy 2-D snapshot files.
if nargin < 6 || isempty(mustExist)
    mustExist = true;
end
if nargin < 5
    targetTime = [];
end
solver = utils.normalize_solver_name(solverName);
Ntheta = max(0, round(double(Ntheta)));
useNx = ~(isempty(Nx) || ~isfinite(double(Nx)));
if useNx
    Nx = max(0, round(double(Nx)));
    patterns = {sprintf('snapshot_%s_Nx%04d_K%04d_t*.mat', solver, Nx, Ntheta), ...
                sprintf('snapshot_%s_*K%04d_t*.mat', solver, Ntheta)};
else
    patterns = {sprintf('snapshot_%s_K%04d_t*.mat', solver, Ntheta), ...
                sprintf('snapshot_%s_*K%04d_t*.mat', solver, Ntheta)};
end
% 新版快照直接在 dataDir 下；同时兼容旧版 dataDir/snapshots。
searchDirs = {char(dataDir), fullfile(char(dataDir), 'snapshots')};
candidates = struct('folder', {}, 'name', {}, 'time', {});
for d = 1:numel(searchDirs)
    if exist(searchDirs{d}, 'dir') ~= 7
        continue;
    end
    for p = 1:numel(patterns)
        files = dir(fullfile(searchDirs{d}, patterns{p}));
        for i = 1:numel(files)
            tt = local_parse_time(files(i).name);
            if isfinite(tt)
                candidates(end+1).folder = files(i).folder; %#ok<AGROW>
                candidates(end).name = files(i).name;
                candidates(end).time = tt;
            end
        end
    end
end
if isempty(candidates)
    fileName = '';
    selectedTime = NaN;
    found = false;
    if mustExist
        if useNx
            error('utils.find_legacy_snapshot_file:NotFound', ...
                '没有找到 legacy snapshot：solver=%s, Nx=%d, K=%d, dataDir=%s', solver, Nx, Ntheta, dataDir);
        else
            error('utils.find_legacy_snapshot_file:NotFound', ...
                '没有找到 legacy snapshot：solver=%s, K=%d, dataDir=%s', solver, Ntheta, dataDir);
        end
    end
    return;
end
fullNames = arrayfun(@(s) fullfile(s.folder, s.name), candidates, 'UniformOutput', false);
[~, ia] = unique(fullNames, 'stable');
candidates = candidates(ia);
times = [candidates.time];
if isempty(targetTime) || ~isfinite(targetTime)
    [~, id] = max(times);
else
    [~, id] = min(abs(times - targetTime));
end
selectedTime = times(id);
fileName = fullfile(candidates(id).folder, candidates(id).name);
found = true;
end

function t = local_parse_time(name)
t = NaN;
tok = regexp(name, '_t([0-9eE+\-.]+)\.mat$', 'tokens', 'once');
if ~isempty(tok)
    t = str2double(tok{1});
end
end
