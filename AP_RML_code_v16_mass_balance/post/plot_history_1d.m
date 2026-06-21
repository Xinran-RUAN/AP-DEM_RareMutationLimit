function figInfo = plot_history_1d(dataDirs, opt)
%PLOT_HISTORY_1D 画一个或多个数据目录的 history 诊断。

if ischar(dataDirs) || local_isstring(dataDirs)
    dataDirs = cellstr(dataDirs);
end
if nargin < 2 || isempty(opt)
    opt = struct();
end
opt = local_fill_defaults(opt);
utils.ensure_dir(opt.outputDir);

fig = figure('Name', 'history diagnostics');
set(fig, 'Color', 'w');

ax1 = subplot(3,1,1); hold(ax1, 'on'); grid(ax1, 'on');
ylabel(ax1, 'residual'); set(ax1, 'YScale', 'log');

ax2 = subplot(3,1,2); hold(ax2, 'on'); grid(ax2, 'on');
ylabel(ax2, '\theta peak');

ax3 = subplot(3,1,3); hold(ax3, 'on'); grid(ax3, 'on');
ylabel(ax3, 'rho max'); xlabel(ax3, 't');

legendText = {};
usedFiles = {};
for i = 1:numel(dataDirs)
    dataDir = dataDirs{i};
    solverName = local_infer_solver(dataDir);
    [dataFile, found] = local_find_history_or_final(dataDir, solverName);
    if ~found
        warning('plot_history_1d:NoData', '跳过目录，没有找到 history/final：%s', dataDir);
        continue;
    end

    S = load(dataFile);
    if ~isfield(S, 'history') || isempty(S.history) || ~isfield(S.history, 't') || isempty(S.history.t)
        warning('plot_history_1d:NoHistory', '文件没有 history：%s', dataFile);
        continue;
    end
    h = S.history;
    t = h.t(:);

    if isfield(h, 'residual')
        semilogy(ax1, t, max(h.residual(:), realmin));
    end

    thetaSeries = [];
    if strcmp(solverName, 'direct') && isfield(h, 'theta_direct')
        thetaSeries = h.theta_direct(:);
    elseif isfield(h, 'theta_wkb')
        thetaSeries = h.theta_wkb(:);
    elseif isfield(h, 'theta_direct')
        thetaSeries = h.theta_direct(:);
    end
    if ~isempty(thetaSeries)
        plot(ax2, t, thetaSeries);
    end

    if isfield(h, 'rho_max')
        plot(ax3, t, h.rho_max(:));
    end

    legendText{end+1} = local_short_label(dataDir); %#ok<AGROW>
    usedFiles{end+1} = dataFile; %#ok<AGROW>
end

if ~isempty(legendText)
    legend(ax1, legendText, 'Interpreter', 'none', 'Location', 'best');
    legend(ax2, legendText, 'Interpreter', 'none', 'Location', 'best');
    legend(ax3, legendText, 'Interpreter', 'none', 'Location', 'best');
end

local_sgtitle('history diagnostics');
figFile = local_save_figure(fig, opt, 'history');

figInfo = struct('job', 'history', 'dataDir', strjoin(dataDirs, ';'), ...
    'dataFile', strjoin(usedFiles, ';'), 'figureFile', figFile);

if logical(opt.closeFigure)
    close(fig);
end
end

function opt = local_fill_defaults(opt)
if ~isfield(opt, 'outputDir') || isempty(opt.outputDir), opt.outputDir = fullfile(utils.project_root(), 'data', 'figures'); end
if ~isfield(opt, 'figurePrefix') || isempty(opt.figurePrefix), opt.figurePrefix = 'post'; end
if ~isfield(opt, 'savePng') || isempty(opt.savePng), opt.savePng = true; end
if ~isfield(opt, 'saveFig') || isempty(opt.saveFig), opt.saveFig = false; end
if ~isfield(opt, 'closeFigure') || isempty(opt.closeFigure), opt.closeFigure = false; end
end

function [dataFile, found] = local_find_history_or_final(dataDir, solverName)
histFile = utils.solver_file(dataDir, solverName, 'history');
if exist(histFile, 'file') == 2
    dataFile = histFile;
    found = true;
    return;
end
[dataFile, found] = utils.find_final_file(dataDir, solverName, false);
end

function solverName = local_infer_solver(dataDir)
[~, name] = fileparts(dataDir);
if startsWith(lower(name), 'direct_') || contains(lower(name), 'direct')
    solverName = 'direct';
else
    solverName = 'wkb';
end
end

function label = local_short_label(dataDir)
[~, label] = fileparts(dataDir);
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

function local_sgtitle(txt)
if exist('sgtitle', 'file') == 2 || exist('sgtitle', 'builtin')
    sgtitle(txt, 'Interpreter', 'none');
else
    annotation('textbox', [0 0.95 1 0.04], 'String', txt, ...
        'Interpreter', 'none', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
end
