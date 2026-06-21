function [S, selectedTime, sourceFile, found] = load_compare_snapshot(dataDir, targetTime, mustExist)
%LOAD_COMPARE_SNAPSHOT Load a paired Test-3 comparison snapshot.
%  Preferred individual files:
%      snapshots_compare/snapshot_compare_###_tT.mat
%  Collection files:
%      snapshots_compare.mat, variable snapshots.
if nargin < 3 || isempty(mustExist)
    mustExist = true;
end
if nargin < 2
    targetTime = [];
end
searchDirs = {fullfile(dataDir, 'snapshots_compare'), dataDir};
patterns = {'snapshot_compare_*_t*.mat'};
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
if ~isempty(candidates)
    fullNames = arrayfun(@(s) fullfile(s.folder, s.name), candidates, 'UniformOutput', false);
    [~, ia] = unique(fullNames, 'stable');
    candidates = candidates(ia);
    times = [candidates.time];
    if isempty(targetTime) || ~isfinite(targetTime)
        [~, id] = max(times);
    else
        [~, id] = min(abs(times - targetTime));
    end
    sourceFile = fullfile(candidates(id).folder, candidates(id).name);
    selectedTime = candidates(id).time;
    L = load(sourceFile);
    S = local_select_snapshot_struct(L, selectedTime);
    found = true;
    return;
end
collectionNames = {char(utils.compare_collection_file(dataDir))};
for i = 1:numel(collectionNames)
    f = collectionNames{i};
    if exist(f, 'file') ~= 2
        continue;
    end
    L = load(f);
    if ~isfield(L, 'snapshots') || isempty(L.snapshots)
        continue;
    end
    times = [L.snapshots.t];
    if isempty(targetTime) || ~isfinite(targetTime)
        [~, id] = max(times);
    else
        [~, id] = min(abs(times - targetTime));
    end
    S = L.snapshots(id);
    selectedTime = times(id);
    sourceFile = f;
    found = true;
    return;
end
S = struct();
selectedTime = NaN;
sourceFile = '';
found = false;
if mustExist
    error('utils.load_compare_snapshot:NotFound', '没有找到 Test-3 compare snapshot。目录：%s', dataDir);
end
end

function S = local_select_snapshot_struct(L, selectedTime)
if isfield(L, 'one')
    S = L.one;
elseif isfield(L, 'snapshots') && ~isempty(L.snapshots)
    times = [L.snapshots.t];
    [~, id] = min(abs(times - selectedTime));
    S = L.snapshots(id);
else
    S = L;
end
end

function t = local_parse_time(name)
t = NaN;
tok = regexp(name, '_t([0-9eE+\-.]+)\.mat$', 'tokens', 'once');
if ~isempty(tok)
    t = str2double(tok{1});
end
end
