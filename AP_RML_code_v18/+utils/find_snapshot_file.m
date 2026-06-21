function [fileName, selectedTime, found] = find_snapshot_file(dataDir, solverName, targetTime, mustExist)
%FIND_SNAPSHOT_FILE Locate the snapshot closest to targetTime.
%
%  新版文件名：
%      <dataDir>/snapshot_<solver>_tT.mat
%
%  同时兼容旧文件名：
%      <dataDir>/snapshot_<solver>_0000123_tT.mat
%      <dataDir>/snapshot_<solver>_step0000123_tT.mat
if nargin < 4 || isempty(mustExist)
    mustExist = true;
end
if nargin < 3
    targetTime = [];
end

solver = utils.normalize_solver_name(solverName);

% 新版快照直接在 dataDir 下；同时兼容旧版 dataDir/snapshots。
searchDirs = {dataDir, fullfile(dataDir, 'snapshots')};

% 优先找新版只带时间的文件，再兼容旧版带编号或 step 的文件。
patterns = {sprintf('snapshot_%s_t*.mat', solver), ...
            sprintf('snapshot_%s_[0-9]*_t*.mat', solver), ...
            sprintf('snapshot_%s_step*_t*.mat', solver)};

candidates = struct('folder', {}, 'name', {}, 'time', {}, 'priority', {});
for d = 1:numel(searchDirs)
    if exist(searchDirs{d}, 'dir') ~= 7
        continue;
    end
    for p = 1:numel(patterns)
        files = dir(fullfile(searchDirs{d}, patterns{p}));
        for i = 1:numel(files)
            tt = local_parse_time(files(i).name);
            if ~isfinite(tt)
                continue;
            end
            candidates(end+1).folder = files(i).folder; %#ok<AGROW>
            candidates(end).name = files(i).name;
            candidates(end).time = tt;
            candidates(end).priority = p;
        end
    end
end

if isempty(candidates)
    fileName = '';
    selectedTime = NaN;
    found = false;
    if mustExist
        error('utils.find_snapshot_file:NotFound', ...
            '没有找到 %s snapshot。期望目录：%s 或 %s', solver, searchDirs{1}, searchDirs{2});
    end
    return;
end

fullNames = arrayfun(@(s) fullfile(s.folder, s.name), candidates, 'UniformOutput', false);
[~, ia] = unique(fullNames, 'stable');
candidates = candidates(ia);

times = [candidates.time];
priorities = [candidates.priority];

if isempty(targetTime) || ~isfinite(targetTime)
    maxTime = max(times);
    ids = find(abs(times - maxTime) <= 10 * eps(max(1, abs(maxTime))));
else
    minErr = min(abs(times - targetTime));
    ids = find(abs(abs(times - targetTime) - minErr) <= 10 * eps(max(1, abs(targetTime))));
end

% 若多个文件对应同一时间，优先使用新版 t-only 文件。
[~, j] = min(priorities(ids));
id = ids(j);

selectedTime = times(id);
fileName = fullfile(candidates(id).folder, candidates(id).name);
found = true;
end

function t = local_parse_time(name)
t = NaN;

% 兼容两类时间片段：
%   旧版：t2.5.mat, t1e-4.mat
%   新版：t2p5.mat, t1em4.mat
% 注意这里不能用 [^.] 截断小数点，否则 t2.5 会被误读。
tok = regexp(name, '_t(.+)\.mat$', 'tokens', 'once');
if isempty(tok)
    return;
end

str = tok{1};
str = strrep(str, 'ep', 'e');
str = strrep(str, 'em', 'e-');
str = strrep(str, 'p', '.');

% 如果是普通负数 m2p5，也做兼容。
if numel(str) > 1 && strncmp(str, 'm', 1) && isempty(strfind(str, 'e'))
    str = ['-' str(2:end)];
end

t = str2double(str);
end
