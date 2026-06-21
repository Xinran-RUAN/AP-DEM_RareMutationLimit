function outdir = find_case_dir(parentDir, preferredCase, pattern, mustExist)
%FIND_CASE_DIR Find a case subdirectory under a data output directory.
%  First tries fullfile(parentDir, preferredCase). If it does not exist, falls
%  back to the newest directory matching pattern, e.g. 'eps_*'.
if nargin < 4 || isempty(mustExist)
    mustExist = true;
end
if nargin < 3 || isempty(pattern)
    pattern = '*';
end
parentDir = char(parentDir);
if nargin >= 2 && ~isempty(preferredCase)
    candidate = fullfile(parentDir, char(preferredCase));
    if exist(candidate, 'dir') == 7
        outdir = candidate;
        return;
    end
end
d = dir(fullfile(parentDir, char(pattern)));
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));
if ~isempty(d)
    [~, id] = max([d.datenum]);
    outdir = fullfile(d(id).folder, d(id).name);
    return;
end
outdir = fullfile(parentDir, char(preferredCase));
if mustExist
    error('utils.find_case_dir:NotFound', ...
        '没有找到 case 数据目录。parent=%s, preferred=%s, pattern=%s', ...
        parentDir, char(preferredCase), char(pattern));
end
end
