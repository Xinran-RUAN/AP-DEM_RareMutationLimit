function figInfo = plot_snapshot_1d(dataDirs, fileNames, opt)
%PLOT_SNAPSHOT_1D 画 final 或指定 snapshot 的 rho、trait marginal 和 n heatmap。

if ischar(dataDirs) || local_isstring(dataDirs)
    dataDirs = cellstr(dataDirs);
end
if nargin < 2 || isempty(fileNames)
    fileNames = {};
elseif ischar(fileNames) || local_isstring(fileNames)
    fileNames = cellstr(fileNames);
end
if nargin < 3 || isempty(opt)
    opt = struct();
end
opt = local_fill_defaults(opt);
utils.ensure_dir(opt.outputDir);

figInfo = struct('job', {}, 'dataDir', {}, 'dataFile', {}, 'figureFile', {});

for i = 1:numel(dataDirs)
    dataDir = dataDirs{i};
    solverName = local_infer_solver(dataDir);
    localFileNames = local_select_file_names(fileNames, i, numel(dataDirs));

    if isempty(localFileNames)
        [dataFile, selectedTime, found] = local_find_snapshot_or_final(dataDir, solverName, opt.snapshotTime);
        if ~found
            warning('plot_snapshot_1d:NoData', '跳过目录，没有找到 final/snapshot：%s', dataDir);
            continue;
        end
        localFileNames = {dataFile};
    else
        selectedTime = NaN;
    end

    for jf = 1:numel(localFileNames)
        dataFile = local_resolve_file(dataDir, localFileNames{jf});
        if exist(dataFile, 'file') ~= 2
            warning('plot_snapshot_1d:MissingFile', '跳过不存在文件：%s', dataFile);
            continue;
        end

        S = load(dataFile);
        data = local_unpack_solution(S, solverName);

        fig = figure('Name', ['snapshot ' local_short_label(dataDir)]);
        set(fig, 'Color', 'w');

        subplot(2,2,1);
        plot(data.x, data.rho, 'LineWidth', 1.2);
        grid on;
        xlabel('x'); ylabel('\rho(x)');
        title(sprintf('%s rho, t=%.6g', solverName, data.t), 'Interpreter', 'none');

        subplot(2,2,2);
        plot(data.theta, data.P, 'LineWidth', 1.2);
        grid on;
        xlabel('\theta'); ylabel('P(\theta)');
        title('trait marginal');

        subplot(2,2,[3 4]);
        imagesc(data.theta, data.x, data.n);
        set(gca, 'YDir', 'normal');
        xlabel('\theta'); ylabel('x');
        title('n(x,\theta)');
        colorbar;

        local_sgtitle(local_title(dataDir, dataFile, data.t, selectedTime));

        suffix = sprintf('%s_%s_t%s', local_safe_token(local_short_label(dataDir)), solverName, local_num_token(data.t));
        figFile = local_save_figure(fig, opt, suffix);

        figInfo(end+1) = struct('job', 'snapshot', 'dataDir', dataDir, ...
            'dataFile', dataFile, 'figureFile', figFile); %#ok<AGROW>

        if logical(opt.closeFigure)
            close(fig);
        end
    end
end
end

function opt = local_fill_defaults(opt)
if ~isfield(opt, 'snapshotTime'), opt.snapshotTime = []; end
if ~isfield(opt, 'outputDir') || isempty(opt.outputDir), opt.outputDir = fullfile(utils.project_root(), 'data', 'figures'); end
if ~isfield(opt, 'figurePrefix') || isempty(opt.figurePrefix), opt.figurePrefix = 'post'; end
if ~isfield(opt, 'savePng') || isempty(opt.savePng), opt.savePng = true; end
if ~isfield(opt, 'saveFig') || isempty(opt.saveFig), opt.saveFig = false; end
if ~isfield(opt, 'closeFigure') || isempty(opt.closeFigure), opt.closeFigure = false; end
end

function names = local_select_file_names(fileNames, index, numDataDirs)
if isempty(fileNames)
    names = {};
    return;
end

% 一个数据目录时，fileNames 可以直接写多个文件名。
if numDataDirs == 1 && index == 1 && numel(fileNames) > 1
    names = fileNames;
    return;
end

if numel(fileNames) == 1
    names = fileNames;
elseif index <= numel(fileNames)
    names = fileNames(index);
else
    names = {};
end

% 支持 fileNames = {{'a.mat','b.mat'}, {'c.mat'}} 这种写法。
if numel(names) == 1 && iscell(names{1})
    names = names{1};
end
end

function [dataFile, selectedTime, found] = local_find_snapshot_or_final(dataDir, solverName, targetTime)
if ~isempty(targetTime) && isnumeric(targetTime) && isscalar(targetTime) && isfinite(targetTime)
    [dataFile, selectedTime, found] = utils.find_snapshot_file(dataDir, solverName, targetTime, false);
    if found
        return;
    end
end
selectedTime = NaN;
[dataFile, found] = utils.find_final_file(dataDir, solverName, false);
end

function dataFile = local_resolve_file(dataDir, fileName)
fileName = char(fileName);
if utils.is_absolute_path(fileName)
    dataFile = fileName;
else
    dataFile = fullfile(dataDir, fileName);
end
end

function solverName = local_infer_solver(dataDir)
[~, name] = fileparts(dataDir);
if startsWith(lower(name), 'direct_') || contains(lower(name), 'direct')
    solverName = 'direct';
else
    solverName = 'wkb';
end
end

function data = local_unpack_solution(S, solverName)
if isfield(S, 'op')
    op = S.op;
else
    op = src.build_operators_1d(S.grid);
end
if isfield(op, 'x')
    x = op.x;
else
    x = S.grid.x;
end
if isfield(op, 'theta')
    theta = op.theta;
else
    theta = S.grid.theta;
end

if isfield(S, 'n') && ~isempty(S.n)
    n = S.n;
elseif isfield(S, 'W') && isfield(S, 'u') && isfield(S, 'par')
    [~, n] = src.reconstruct_rho_1d(S.W, S.u, S.par.eps, op, 'direct-log');
else
    error('plot_snapshot_1d:NoDensity', '文件中既没有 n，也没有 W/u 可重构 n。');
end

if isfield(S, 'rho') && ~isempty(S.rho)
    rho = S.rho;
else
    rho = op.dtheta * sum(n, 2);
end

wx = utils.trapz_weights(op.nx, op.dx);
P = wx.' * n;

if isfield(S, 't')
    t = S.t;
else
    t = NaN;
end

data = struct();
data.x = x;
data.theta = theta;
data.n = n;
data.rho = rho;
data.P = P;
data.t = t;
data.solver = solverName;
end

function titleText = local_title(dataDir, dataFile, t, selectedTime)
[~, dirName] = fileparts(dataDir);
[~, fileBase, ext] = fileparts(dataFile);
if isfinite(selectedTime)
    titleText = sprintf('%s | %s%s | t=%.6g, requested %.6g', dirName, fileBase, ext, t, selectedTime);
else
    titleText = sprintf('%s | %s%s | t=%.6g', dirName, fileBase, ext, t);
end
end

function figFile = local_save_figure(fig, opt, suffix)
base = fullfile(opt.outputDir, [local_safe_token(opt.figurePrefix) '_' suffix]);
figFile = '';
if logical(opt.savePng)
    figFile = [base '.png'];
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, figFile, 'Resolution', 200);
    else
        saveas(fig, figFile);
    end
end
if logical(opt.saveFig)
    savefig(fig, [base '.fig']);
end
end

function label = local_short_label(dataDir)
[~, label] = fileparts(dataDir);
end

function tf = local_isstring(v)
tf = false;
if exist('isstring', 'builtin') || exist('isstring', 'file')
    tf = isstring(v);
end
end

function s = local_to_char(v)
if ischar(v)
    s = v;
elseif local_isstring(v)
    s = char(v);
elseif isnumeric(v) && isscalar(v)
    s = num2str(v);
else
    try
        s = char(v);
    catch
        s = 'figure';
    end
end
end

function s = local_safe_token(v)
s = local_to_char(v);
s = regexprep(s, '[^A-Za-z0-9_\-]+', '_');
s = regexprep(s, '^_+|_+$', '');
if isempty(s), s = 'figure'; end
end

function s = local_num_token(x)
if isempty(x) || ~isnumeric(x) || ~isscalar(x) || ~isfinite(x)
    s = 'nan';
    return;
end
x = double(x);
if abs(x - round(x)) < 1e-12 && abs(x) < 1e8
    s = sprintf('%g', x);
else
    s = sprintf('%.6g', x);
end
s = lower(s);
s = strrep(s, '+', '');
s = regexprep(s, 'e-0*', 'em');
s = regexprep(s, 'e\+?0*', 'ep');
s = strrep(s, '.', 'p');
s = strrep(s, '-', 'm');
s = regexprep(s, '[^A-Za-z0-9]+', '_');
end

function local_sgtitle(txt)
if exist('sgtitle', 'file') == 2 || exist('sgtitle', 'builtin')
    sgtitle(txt, 'Interpreter', 'none');
else
    annotation('textbox', [0 0.95 1 0.04], 'String', txt, ...
        'Interpreter', 'none', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
end
