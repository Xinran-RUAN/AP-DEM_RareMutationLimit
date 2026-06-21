function [fileName, found] = find_checkpoint_file(dataDir, solverName, mustExist)
%FIND_CHECKPOINT_FILE Locate canonical latest checkpoint.
if nargin < 3 || isempty(mustExist)
    mustExist = true;
end
solver = utils.normalize_solver_name(solverName);
fileName = utils.solver_file(dataDir, solver, 'checkpoint');
found = exist(fileName, 'file') == 2;
if ~found
    if mustExist
        error('utils.find_checkpoint_file:NotFound', ...
            '没有找到 %s checkpoint 文件。期望文件：%s', solver, fileName);
    end
    fileName = '';
end
end
