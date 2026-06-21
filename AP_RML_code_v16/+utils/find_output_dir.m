function outdir = find_output_dir(names, mustExist)
%FIND_OUTPUT_DIR Return the first existing data/<name> directory.
%  names may be a char/string or a cell array of aliases ordered by priority.
if nargin < 2 || isempty(mustExist)
    mustExist = true;
end
if nargin < 1 || isempty(names)
    names = {''};
elseif ischar(names) || isstring(names)
    names = cellstr(names);
end
candidates = cell(size(names));
for i = 1:numel(names)
    nm = char(names{i});
    if isempty(nm)
        candidates{i} = utils.output_dir();
    elseif utils.is_absolute_path(nm)
        candidates{i} = nm;
    else
        candidates{i} = utils.output_dir(nm);
    end
    if exist(candidates{i}, 'dir') == 7
        outdir = candidates{i};
        return;
    end
end
outdir = candidates{1};
if mustExist
    msg = sprintf('  %s\n', candidates{:});
    error('utils.find_output_dir:NotFound', '没有找到任何候选数据目录：\n%s', msg);
end
end
