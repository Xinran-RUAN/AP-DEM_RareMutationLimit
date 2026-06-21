function figInfo = plot_case_1d(dataDirs, fileNames, opt)
%PLOT_CASE_1D 一维数据后处理总入口。
%
% 输入：
%   dataDirs  : 一个或多个数据目录。
%   fileNames : 可为空；也可以给一个或多个 mat 文件名。
%   opt.jobs  : {'history','snapshot'}。
%
% 示例：
%   plot_case_1d({dataDir1,dataDir2}, {}, struct('jobs',{{'history','snapshot'}}));

if nargin < 1 || isempty(dataDirs)
    error('plot_case_1d:NoDataDirs', '请提供一个或多个数据目录。');
end
if nargin < 2
    fileNames = {};
end
if nargin < 3 || isempty(opt)
    opt = struct();
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

if ischar(dataDirs) || local_isstring(dataDirs)
    dataDirs = cellstr(dataDirs);
end
if ischar(fileNames) || local_isstring(fileNames)
    fileNames = cellstr(fileNames);
end

opt = local_fill_plot_defaults(opt, root);
utils.ensure_dir(opt.outputDir);

figInfo = struct('job', {}, 'dataDir', {}, 'dataFile', {}, 'figureFile', {});

jobs = opt.jobs;
if ischar(jobs) || local_isstring(jobs)
    jobs = cellstr(jobs);
end

if any(strcmpi(jobs, 'history'))
    hInfo = plot_history_1d(dataDirs, opt);
    figInfo = [figInfo, hInfo]; %#ok<AGROW>
end

if any(strcmpi(jobs, 'snapshot')) || any(strcmpi(jobs, 'final'))
    sInfo = plot_snapshot_1d(dataDirs, fileNames, opt);
    figInfo = [figInfo, sInfo]; %#ok<AGROW>
end
end

function opt = local_fill_plot_defaults(opt, root)
if ~isfield(opt, 'jobs') || isempty(opt.jobs), opt.jobs = {'history','snapshot'}; end
if ~isfield(opt, 'snapshotTime'), opt.snapshotTime = []; end
if ~isfield(opt, 'outputDir') || isempty(opt.outputDir), opt.outputDir = fullfile(root, 'data', 'figures'); end
if ~isfield(opt, 'figurePrefix') || isempty(opt.figurePrefix), opt.figurePrefix = 'post'; end
if ~isfield(opt, 'savePng') || isempty(opt.savePng), opt.savePng = true; end
if ~isfield(opt, 'saveFig') || isempty(opt.saveFig), opt.saveFig = false; end
if ~isfield(opt, 'closeFigure') || isempty(opt.closeFigure), opt.closeFigure = false; end
end

function tf = local_isstring(v)
tf = false;
if exist('isstring', 'builtin') || exist('isstring', 'file')
    tf = isstring(v);
end
end
