function [fileName, found] = find_final_file(dataDir, solverName, mustExist)
%FIND_FINAL_FILE Locate a final result file, preferring canonical names.
if nargin < 3 || isempty(mustExist)
    mustExist = true;
end
solver = utils.normalize_solver_name(solverName);
canon = utils.solver_file(dataDir, solver, 'final');
if exist(canon, 'file') == 2
    fileName = canon;
    found = true;
    return;
end
switch solver
    case 'direct'
        patterns = {'result_direct_final.mat', '*_direct_final.mat'};
    case 'wkb'
        patterns = {'result_wkb_final.mat', '*_final.mat'};
    otherwise
        patterns = {sprintf('result_%s_final.mat', solver), '*_final.mat'};
end
files = [];
for p = 1:numel(patterns)
    f = dir(fullfile(dataDir, patterns{p}));
    if strcmp(solver, 'wkb')
        f = f(~contains({f.name}, '_direct_final'));
        f = f(~contains({f.name}, 'result_direct_final'));
    end
    files = [files; f(:)]; %#ok<AGROW>
end
if isempty(files)
    fileName = '';
    found = false;
    if mustExist
        error('utils.find_final_file:NotFound', '没有找到 %s final 文件。目录：%s', solver, dataDir);
    end
    return;
end
[~, id] = max([files.datenum]);
fileName = fullfile(files(id).folder, files(id).name);
found = true;
end
