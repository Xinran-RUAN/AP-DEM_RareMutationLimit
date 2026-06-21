%% POST_MAIN_1D
% -------------------------------------------------------------------------
% 一维统一后处理：
%   自动识别 direct / WKB 数据，并采用相应方式构造 P(theta)。
%   默认扫描 data/caseName 下的 result_*_final.mat。
%
% direct:
%   P_k = int_x n(x,theta_k) dx,
%   log P_k = log(P_k).
%
% WKB:
%   A_k = int_x W(x,theta_k) dx,
%   log P_k = u_k/eps + log(A_k).
%
% 默认推荐：
%   用 WENO 只定位峰值位置 theta_*，
%   再固定 theta_* 对 log P 做约束二次拟合，
%   最后用原始粗网格总质量定标。
% -------------------------------------------------------------------------

clear; clc; close all;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

caseName = 'main_1d';
caseDir = fullfile(root, 'data', caseName);

% -------------------------------------------------------------------------
% 数据选择
% -------------------------------------------------------------------------
% 推荐默认用 auto：扫描 data/caseName 的所有子目录，自动寻找
% result_wkb_final.mat / result_direct_final.mat，并根据文件名或变量字段识别类型。
%
% 若只想画指定数据，把 post.dataSelectMode 改为 'manual'，并在 files 中写条目。
% 每个条目可写：
%   struct('dirName','...', 'fileName','auto', 'targetTime',10, 'solver','auto', 'label','...')
%
% fileName = 'auto' 时：
%   targetTime 为空  -> 优先找 final result；
%   targetTime 非空 -> 找最接近 targetTime 的 snapshot。
% -------------------------------------------------------------------------

post.dataSelectMode = 'auto';     % {'auto','manual'}
post.autoFileMode   = 'final';    % {'final','snapshot'}; auto 模式使用
post.autoTargetTime = [];         % autoFileMode='snapshot' 时填写，例如 10
post.autoSolvers    = {'wkb','direct'};   % {'wkb'}, {'direct'}, {'wkb','direct'}
post.autoMaxItems   = Inf;        % 自动扫描最多画多少条曲线
post.autoSkipFigureDir = true;    % 自动扫描时跳过 data/figures

% manual 模式示例。auto 模式下可保持 files = {}。
files = { ...
       'fileName','snapshot_wkb_t10.mat', ...
       'solver','wkb', ...
       'label','WKB, N_W=16, N_u=64', ...   
};
% files = { ...
%     struct('dirName','wkb_eps1epm5_dt0p001_Nx50_Ntheta16_t10_Nu64', ...
%            'fileName','result_wkb_final.mat', ...
%            'solver','wkb', ...
%            'label','WKB, N_W=16, N_u=64'), ...
%     struct('dirName','direct_eps1epm5_dt0p001_Nx50_Ntheta64_t10', ...
%            'fileName','result_direct_final.mat', ...
%            'solver','direct', ...
%            'label','Direct, N_\theta=64') ...
%     };

if strcmpi(post.dataSelectMode, 'auto') || isempty(files)
    files = local_collect_auto_file_items(caseDir, post);
end

% -------------------------------------------------------------------------
% P(theta) 重构方式
% -------------------------------------------------------------------------

post.PReconMode = 'weno-anchored-quadratic-mass';
% {'weno-anchored-quadratic-mass','log-quadratic-mass','log-pchip','pointwise-exp'}

post.PMassMode = 'raw-data';
% {'raw-data','none','unit'}
%   raw-data:
%       默认推荐。用总质量定标，不是 normalize。若 useSolverCorrectedMassIfAvailable=true
%       且文件中含 massInfo，则优先使用求解过程质量修正后的质量；否则使用粗网格 P 质量。
%   none:
%       不做质量定标，直接使用拟合得到的绝对 logP 高度。
%   unit:
%       归一化为 int P dtheta = 1，仅用于比较形状。

% 如果求解过程中启用了质量修正，snapshot/result 中会保存 massInfo。
% 这里为 true 时，后处理会优先使用求解器修正后的质量作为可信总质量；
% 若文件中没有 massInfo，则自动退回 PMassMode 指定的 raw-data/unit/none。
post.useSolverCorrectedMassIfAvailable = true;   % {true,false}

post.NthetaPlot = 4096;

% -------------------------------------------------------------------------
% 用于峰值定位
% -------------------------------------------------------------------------

post.peakLocationMethod = 'weno5';   % {'weno5','pchip','spline','linear','makima','raw'}
post.peakSearchNtheta = 4096;
post.peakRefineParabolic = true;

% -------------------------------------------------------------------------
% 用于峰附近 log P 拟合
% -------------------------------------------------------------------------

post.peakFitPoints = 5;
post.peakFitDegree = 2;

post.peakFitWeighted = true;
post.peakFitWeightMode = 'log';      % {'none','log','gaussian'}
post.peakFitLogDrop = 8;

post.logInterpMethod = 'pchip';      % {'pchip','spline','linear','makima'}

% pointwise-exp 模式才会用到
post.directThetaInterp = 'pchip';    % {'weno5','pchip','spline','linear','makima'}
post.wkbThetaInterp    = 'pchip';    % {'weno5','pchip','spline','linear','makima'}
post.clipNegative = true;

% -------------------------------------------------------------------------
% 画图选项
% -------------------------------------------------------------------------

post.showRawPoints = true;
post.maxRawMarkers = 8;
post.rawMarkerMode = 'uniform_peak'; % {'uniform','uniform_peak','all'}

post.lineWidth = 2.4;
post.rawMarkerSize = 6.5;
post.rawMarkerLineWidth = 1.3;

post.plotLogDiagnostic = true;
post.peakWindowPlot = true;
post.plotPeakLocation = true;

post.outputDir = fullfile(root, 'data', 'figures', caseName);
post.figurePrefix = 'post_main_1d_Ptheta';

if exist(post.outputDir, 'dir') ~= 7
    mkdir(post.outputDir);
end

%% ========================= 计算并画 P(theta) =========================

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [4, 4, 18, 12]);

ax = axes(fig);
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

fprintf('\n================ log P(theta) 后处理 ================\n');
fprintf('PReconMode          = %s\n', post.PReconMode);
fprintf('PMassMode           = %s\n', post.PMassMode);
fprintf('use solver mass     = %d\n', post.useSolverCorrectedMassIfAvailable);
fprintf('peakLocationMethod  = %s\n', post.peakLocationMethod);
fprintf('NthetaPlot          = %d\n', post.NthetaPlot);
fprintf('peakSearchNtheta    = %d\n', post.peakSearchNtheta);
fprintf('peakFitPoints       = %d\n', post.peakFitPoints);
fprintf('data items          = %d\n', numel(files));

allData = cell(numel(files), 1);

for q = 1:numel(files)

    item = local_normalize_file_item(files{q}, caseDir);
    [filePath, solverFromFile] = local_resolve_file_path(item);

    Sraw = load(filePath);
    S = local_select_data_struct(Sraw);

    solver = local_infer_solver(item, solverFromFile, S, filePath);

    one = local_compute_Ptheta_log_reconstruction(S, solver, post);
    one.label = item.label;
    one.filePath = filePath;
    one.solver = solver;

    allData{q} = one;

    [lineStyle, markerStyle] = local_plot_style(q, item);

    hLine = plot(ax, one.thetaPlot, one.PPlot, ...
        'LineStyle', lineStyle, ...
        'Marker', 'none', ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', item.label);

    if post.showRawPoints
        markerIndex = local_select_raw_markers( ...
            one.thetaRaw, one.PRaw, post.maxRawMarkers, post.rawMarkerMode);

        if ~isempty(markerIndex)
            plot(ax, one.thetaRaw(markerIndex), one.PRaw(markerIndex), ...
                'LineStyle', 'none', ...
                'Marker', markerStyle, ...
                'MarkerSize', post.rawMarkerSize, ...
                'LineWidth', post.rawMarkerLineWidth, ...
                'Color', hLine.Color, ...
                'MarkerFaceColor', 'none', ...
                'HandleVisibility', 'off');
        end
    end

    fprintf('\n[%s]\n', item.label);
    fprintf('  file            = %s\n', filePath);
    fprintf('  solver          = %s\n', solver);
    fprintf('  eps             = %.16g\n', one.epsVal);
    fprintf('  Nx              = %d\n', one.Nx);
    fprintf('  NthetaRaw       = %d\n', one.NthetaRaw);
    fprintf('  raw mass        = %.16e\n', one.massRaw);
    if isfield(one, 'massForScaling')
        fprintf('  scale mass      = %.16e  (%s)\n', one.massForScaling, one.massSource);
    end
    fprintf('  rec mass        = %.16e\n', one.massPlot);
    fprintf('  raw peak        = theta %.8f, P %.8e\n', ...
        one.thetaPeakRaw, one.maxPRaw);
    fprintf('  rec peak        = theta %.8f, P %.8e\n', ...
        one.thetaPeakPlot, one.maxPPlot);

    if isfield(one, 'fitThetaPeak')
        fprintf('  anchored peak   = %.8f\n', one.fitThetaPeak);
        fprintf('  fit center      = %.8f\n', one.fitThetaCenter);
        fprintf('  fit curvature   = %.8e\n', one.fitCurvature);
        fprintf('  fit fallback    = %s\n', one.fitFallback);
    end
end

xlabel(ax, '\theta', 'Interpreter', 'tex');
ylabel(ax, 'P(\theta)', 'Interpreter', 'tex');
title(ax, 'Trait marginal reconstruction by log P', 'Interpreter', 'tex');
legend(ax, 'Location', 'best', 'Box', 'off');
xlim(ax, [0, 1]);
set(ax, 'FontSize', 12);

pngFile = fullfile(post.outputDir, [post.figurePrefix, '_Ptheta.png']);
figFile = fullfile(post.outputDir, [post.figurePrefix, '_Ptheta.fig']);

try
    exportgraphics(fig, pngFile, 'Resolution', 300);
catch
    saveas(fig, pngFile);
end
savefig(fig, figFile);

fprintf('\nP(theta) 图已保存：\n%s\n%s\n', pngFile, figFile);

%% ========================= log P 诊断图 =========================

if post.plotLogDiagnostic

    figLog = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [5, 5, 18, 12]);

    axLog = axes(figLog);
    hold(axLog, 'on');
    box(axLog, 'on');
    grid(axLog, 'on');

    for q = 1:numel(allData)

        one = allData{q};
        item = local_normalize_file_item(files{q}, caseDir);
        [lineStyle, markerStyle] = local_plot_style(q, item);

        hLine = plot(axLog, one.thetaPlot, one.logPPlot, ...
            'LineStyle', lineStyle, ...
            'Marker', 'none', ...
            'LineWidth', post.lineWidth, ...
            'DisplayName', one.label);

        markerIndex = local_select_raw_markers( ...
            one.thetaRaw, one.logPRaw, post.maxRawMarkers, post.rawMarkerMode);

        if ~isempty(markerIndex)
            plot(axLog, one.thetaRaw(markerIndex), one.logPRaw(markerIndex), ...
                'LineStyle', 'none', ...
                'Marker', markerStyle, ...
                'MarkerSize', post.rawMarkerSize, ...
                'LineWidth', post.rawMarkerLineWidth, ...
                'Color', hLine.Color, ...
                'MarkerFaceColor', 'none', ...
                'HandleVisibility', 'off');
        end

        if post.peakWindowPlot && isfield(one, 'fitIndex')
            plot(axLog, one.thetaRaw(one.fitIndex), one.logPRaw(one.fitIndex), ...
                'LineStyle', 'none', ...
                'Marker', markerStyle, ...
                'MarkerSize', post.rawMarkerSize + 2, ...
                'LineWidth', post.rawMarkerLineWidth + 0.6, ...
                'Color', hLine.Color, ...
                'MarkerFaceColor', hLine.Color, ...
                'HandleVisibility', 'off');
        end

        if post.plotPeakLocation && isfield(one, 'fitThetaPeak')
            yl = ylim(axLog);
            plot(axLog, [one.fitThetaPeak, one.fitThetaPeak], yl, ...
                'LineStyle', ':', ...
                'LineWidth', 1.2, ...
                'Color', hLine.Color, ...
                'HandleVisibility', 'off');
        end
    end

    xlabel(axLog, '\theta', 'Interpreter', 'tex');
    ylabel(axLog, 'log P(\theta)', 'Interpreter', 'tex');
    title(axLog, 'Diagnostic plot for log P reconstruction', 'Interpreter', 'tex');
    legend(axLog, 'Location', 'best', 'Box', 'off');
    xlim(axLog, [0, 1]);
    set(axLog, 'FontSize', 12);

    pngFile = fullfile(post.outputDir, [post.figurePrefix, '_logP.png']);
    figFile = fullfile(post.outputDir, [post.figurePrefix, '_logP.fig']);

    try
        exportgraphics(figLog, pngFile, 'Resolution', 300);
    catch
        saveas(figLog, pngFile);
    end
    savefig(figLog, figFile);

    fprintf('\nlogP 诊断图已保存：\n%s\n%s\n', pngFile, figFile);
end

%% ========================================================================
% 局部函数
% ========================================================================

function files = local_collect_auto_file_items(caseDir, post)

autoSolvers = local_normalize_solver_list(local_get_post_field(post, 'autoSolvers', {'wkb','direct'}));
autoFileMode = lower(local_get_post_field(post, 'autoFileMode', 'final'));
autoTargetTime = local_get_post_field(post, 'autoTargetTime', []);
autoMaxItems = local_get_post_field(post, 'autoMaxItems', Inf);
skipFigureDir = local_get_post_field(post, 'autoSkipFigureDir', true);

if exist(caseDir, 'dir') ~= 7
    error('caseDir 不存在：\n%s', caseDir);
end

dirs = dir(caseDir);
dirs = dirs([dirs.isdir]);
items = {};

for i = 1:numel(dirs)
    name = dirs(i).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue;
    end
    if skipFigureDir && strcmpi(name, 'figures')
        continue;
    end

    dataDir = fullfile(dirs(i).folder, name);

    for sId = 1:numel(autoSolvers)
        sol = autoSolvers{sId};
        filePath = '';
        found = false;

        switch autoFileMode
            case 'final'
                [filePath, found] = local_try_find_final_file(dataDir, sol);
            case 'snapshot'
                if isempty(autoTargetTime)
                    error('post.autoFileMode=''snapshot'' 时必须设置 post.autoTargetTime。');
                end
                try
                    [filePath, solFound] = local_find_snapshot_by_time_auto(dataDir, autoTargetTime, sol);
                    found = true;
                    sol = solFound;
                catch
                    found = false;
                end
            otherwise
                error('未知 post.autoFileMode：%s。', autoFileMode);
        end

        if found && exist(filePath, 'file') == 2
            [~, fn, ext] = fileparts(filePath);
            it = struct();
            it.dirName = name;
            it.dataDir = dataDir;
            it.fileName = [fn, ext];
            it.targetTime = [];
            it.solver = sol;
            it.label = local_auto_label(name, sol, [fn, ext]);
            items{end+1} = it; %#ok<AGROW>
        end
    end
end

if isempty(items)
    error(['自动扫描没有找到可后处理的数据。\ncaseDir = %s\n' ...
           '请先运行 run_wkb_main_1d.m 或 run_direct_main_1d.m，' ...
           '或者把 post.dataSelectMode 改为 manual 并手动指定 files。'], caseDir);
end

% 排序：先按目录名，再按 solver，保证图例稳定。
keys = cell(numel(items),1);
for i = 1:numel(items)
    keys{i} = sprintf('%s__%s__%s', items{i}.dirName, items{i}.solver, items{i}.fileName);
end
[~, ord] = sort(keys);
items = items(ord);

if isfinite(autoMaxItems) && numel(items) > autoMaxItems
    items = items(1:autoMaxItems);
end

files = items;
end

function val = local_get_post_field(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
end
end

function solvers = local_normalize_solver_list(raw)
if ischar(raw) || isstring(raw)
    raw = {char(raw)};
end
solvers = cell(size(raw));
for i = 1:numel(raw)
    sol = lower(char(raw{i}));
    if strcmp(sol, 'all')
        solvers = {'wkb','direct'};
        return;
    end
    if ~ismember(sol, {'wkb','direct'})
        error('autoSolvers 只能包含 wkb/direct/all，收到：%s', sol);
    end
    solvers{i} = sol;
end
end

function label = local_auto_label(dirName, solver, fileName)
label = sprintf('%s: %s', upper(solver), dirName);
if contains(lower(fileName), 'snapshot')
    label = sprintf('%s (%s)', label, fileName);
end
label = strrep(label, '_', '\_');
end

function [filePath, solver] = local_find_final_file_auto(dataDir, solverHint)
if nargin < 2 || isempty(solverHint)
    solverHint = 'auto';
end
if strcmpi(solverHint, 'auto')
    [fw, okW] = local_try_find_final_file(dataDir, 'wkb');
    [fd, okD] = local_try_find_final_file(dataDir, 'direct');
    if okW && ~okD
        filePath = fw; solver = 'wkb'; return;
    elseif okD && ~okW
        filePath = fd; solver = 'direct'; return;
    elseif okW && okD
        error('目录同时含 WKB 和 direct final。请显式设置 solver。目录：%s', dataDir);
    else
        error('没有找到 final result 文件。目录：%s', dataDir);
    end
else
    [filePath, ok] = local_try_find_final_file(dataDir, solverHint);
    if ~ok
        error('没有找到 %s final result 文件。目录：%s', solverHint, dataDir);
    end
    solver = lower(solverHint);
end
end

function [filePath, found] = local_try_find_final_file(dataDir, solver)
solver = lower(char(solver));
found = false;
filePath = '';

switch solver
    case 'wkb'
        patterns = {'result_wkb_final.mat', '*wkb*final*.mat'};
        reject = {'direct'};
    case 'direct'
        patterns = {'result_direct_final.mat', '*direct*final*.mat'};
        reject = {};
    otherwise
        error('未知 solver：%s。', solver);
end

cand = [];
for p = 1:numel(patterns)
    f = dir(fullfile(dataDir, patterns{p}));
    if ~isempty(reject)
        names = lower({f.name});
        keep = true(size(f));
        for r = 1:numel(reject)
            keep = keep & ~contains(names, reject{r});
        end
        f = f(keep);
    end
    cand = [cand; f(:)]; %#ok<AGROW>
end

if isempty(cand)
    return;
end

[~, id] = max([cand.datenum]);
filePath = fullfile(cand(id).folder, cand(id).name);
found = true;
end

function item = local_normalize_file_item(rawItem, caseDir)

if ischar(rawItem) || isstring(rawItem)
    p = char(rawItem);

    if contains(p, filesep) || endsWith(p, '.mat')
        [folderName, fileName, ext] = fileparts(p);

        if isempty(ext)
            item.dirName = p;
            item.fileName = 'auto';
        else
            item.dirName = folderName;
            item.fileName = [fileName, ext];
        end
    else
        item.dirName = p;
        item.fileName = 'auto';
    end

    item.targetTime = [];
    item.label = local_label_from_path(p);
    item.solver = 'auto';
else
    item = rawItem;

    if ~isfield(item, 'fileName') || isempty(item.fileName)
        item.fileName = 'auto';
    end

    if ~isfield(item, 'targetTime')
        item.targetTime = [];
    end

    if ~isfield(item, 'solver') || isempty(item.solver)
        item.solver = 'auto';
    end

    if ~isfield(item, 'label') || isempty(item.label)
        item.label = local_label_from_path(item.dirName);
    end
end

if isfield(item, 'dataDir') && ~isempty(item.dataDir)
    item.dataDir = item.dataDir;
elseif isfield(item, 'dirName') && ~isempty(item.dirName)
    if exist(item.dirName, 'dir') == 7
        item.dataDir = item.dirName;
    else
        item.dataDir = fullfile(caseDir, item.dirName);
    end
else
    error('每个数据项必须提供 dirName 或 dataDir。');
end
end

function label = local_label_from_path(p)

[~, name, ext] = fileparts(char(p));

if isempty(name)
    name = char(p);
end

label = [name, ext];
label = strrep(label, '_', '\_');
end

function [filePath, solverFromFile] = local_resolve_file_path(item)

solverFromFile = 'auto';

if strcmpi(item.fileName, 'auto')
    if isfield(item, 'targetTime') && ~isempty(item.targetTime)
        [filePath, solverFromFile] = local_find_snapshot_by_time_auto( ...
            item.dataDir, item.targetTime, item.solver);
    else
        [filePath, solverFromFile] = local_find_final_file_auto( ...
            item.dataDir, item.solver);
    end
else
    filePath = fullfile(item.dataDir, item.fileName);

    if exist(filePath, 'file') ~= 2
        error('找不到文件：\n%s', filePath);
    end

    solverFromFile = local_infer_solver_from_name(filePath);
end
end

function solver = local_infer_solver(item, solverFromFile, S, filePath)

solver = 'auto';

if isfield(item, 'solver') && ~isempty(item.solver) ...
        && ~strcmpi(item.solver, 'auto')
    solver = lower(item.solver);
end

if strcmpi(solver, 'auto')
    solver = local_infer_solver_from_name(filePath);
end

if strcmpi(solver, 'auto') && ~strcmpi(solverFromFile, 'auto')
    solver = solverFromFile;
end

if strcmpi(solver, 'auto')
    solver = local_infer_solver_from_struct(S);
end

if ~ismember(lower(solver), {'direct','wkb'})
    error('无法识别数据类型。请在 files 中显式写 solver=''direct'' 或 solver=''wkb''。\n文件：%s', filePath);
end

solver = lower(solver);
end

function solver = local_infer_solver_from_name(name)

name = lower(char(name));

if contains(name, 'direct')
    solver = 'direct';
elseif contains(name, 'wkb')
    solver = 'wkb';
else
    solver = 'auto';
end
end

function solver = local_infer_solver_from_struct(S)

if isfield(S, 'W') && isfield(S, 'u') ...
        && ~isempty(S.W) && ~isempty(S.u)
    solver = 'wkb';
elseif isfield(S, 'n') && ~isempty(S.n)
    solver = 'direct';
else
    solver = 'auto';
end
end

function [filePath, solver] = local_find_snapshot_by_time_auto(dataDir, targetTime, solverHint)

if nargin < 3 || isempty(solverHint)
    solverHint = 'auto';
end

if exist(dataDir, 'dir') ~= 7
    error('数据目录不存在：\n%s', dataDir);
end

if strcmpi(solverHint, 'direct')
    d = dir(fullfile(dataDir, 'snapshot_direct*.mat'));
elseif strcmpi(solverHint, 'wkb')
    d = dir(fullfile(dataDir, 'snapshot_wkb*.mat'));
else
    d1 = dir(fullfile(dataDir, 'snapshot_direct*.mat'));
    d2 = dir(fullfile(dataDir, 'snapshot_wkb*.mat'));
    d = [d1; d2];
end

if isempty(d)
    error('目录中没有找到 snapshot 文件：\n%s', dataDir);
end

times = nan(numel(d), 1);

for i = 1:numel(d)
    times(i) = local_parse_time_from_snapshot_name(d(i).name);
end

valid = isfinite(times);

if ~any(valid)
    error('找到了 snapshot 文件，但无法从文件名解析时间：\n%s', dataDir);
end

idxValid = find(valid);
[~, k] = min(abs(times(valid) - targetTime));
idx = idxValid(k);

filePath = fullfile(d(idx).folder, d(idx).name);
solver = local_infer_solver_from_name(filePath);
end

function t = local_parse_time_from_snapshot_name(name)

t = NaN;

token = regexp(name, '_t([^_\.]+)\.mat$', 'tokens', 'once');

if isempty(token)
    return;
end

str = token{1};
str = strrep(str, 'p', '.');
str = strrep(str, 'm', '-');

t = str2double(str);
end

function S = local_select_data_struct(L)

S = L;

if isfield(L, 'snapshot') && isstruct(L.snapshot)
    S = L.snapshot;
elseif isfield(L, 'direct') && isstruct(L.direct)
    S = L.direct;
elseif isfield(L, 'wkb') && isstruct(L.wkb)
    S = L.wkb;
elseif isfield(L, 'one') && isstruct(L.one)
    S = L.one;
else
    names = fieldnames(L);
    if numel(names) == 1 && isstruct(L.(names{1}))
        S = L.(names{1});
    end
end
end

function [massForScaling, massSource] = local_mass_for_scaling(S, thetaRaw, logPRaw, post)
% 后处理定标质量。若求解中启用了质量修正，优先使用 snapshot/result 里保存的修正质量。
massRaw = local_periodic_integral_log(thetaRaw, logPRaw);
massForScaling = massRaw;
massSource = 'raw-data';

useSolverMass = false;
if isstruct(post) && isfield(post, 'useSolverCorrectedMassIfAvailable') && ~isempty(post.useSolverCorrectedMassIfAvailable)
    useSolverMass = logical(post.useSolverCorrectedMassIfAvailable);
end
if ~useSolverMass
    return;
end

[mSolver, ok, source] = local_get_solver_corrected_mass(S);
if ok
    massForScaling = mSolver;
    massSource = source;
end
end

function [m, ok, source] = local_get_solver_corrected_mass(S)
m = NaN;
ok = false;
source = '';

if isfield(S, 'massInfo') && isstruct(S.massInfo)
    mi = S.massInfo;
    enabled = local_get_bool_field(mi, 'enabled', false);
    if enabled
        candidates = {'correctedMass','massAfter'};
        for i = 1:numel(candidates)
            f = candidates{i};
            if isfield(mi, f) && isnumeric(mi.(f)) && isscalar(mi.(f)) && isfinite(mi.(f)) && mi.(f) > 0
                m = double(mi.(f));
                ok = true;
                source = ['solver-corrected:massInfo.' f];
                return;
            end
        end
    end
end

if isfield(S, 'diagnostics') && isstruct(S.diagnostics)
    dg = S.diagnostics;
    enabled = local_get_bool_field(dg, 'mass_correction_enabled', false);
    if enabled
        candidates = {'mass_corrected','mass_total'};
        for i = 1:numel(candidates)
            f = candidates{i};
            if isfield(dg, f) && isnumeric(dg.(f)) && isscalar(dg.(f)) && isfinite(dg.(f)) && dg.(f) > 0
                m = double(dg.(f));
                ok = true;
                source = ['solver-corrected:diagnostics.' f];
                return;
            end
        end
    end
end
end

function tf = local_get_bool_field(s, name, defaultValue)
tf = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end

function one = local_compute_Ptheta_log_reconstruction(S, solver, post)

[x, thetaRaw] = local_get_grid(S);

Nx = numel(x);
NthetaRaw = numel(thetaRaw);

thetaPlot = linspace(0, 1, post.NthetaPlot + 1).';
thetaPlot(end) = [];

epsVal = local_get_eps_if_possible(S);
fitInfo = struct();

switch lower(solver)

    case 'direct'
        if ~isfield(S, 'n') || isempty(S.n)
            error('direct 数据中没有 n。');
        end

        nRaw = local_to_xt_matrix(S.n, Nx, NthetaRaw, 'n');

        PRaw = trapz(x(:), nRaw, 1).';
        PRaw = real(PRaw(:));
        PRaw(PRaw < 0) = 0;

        logPRaw = log(max(PRaw, realmin));
        [massForScaling, massSource] = local_mass_for_scaling(S, thetaRaw, logPRaw, post);

        if strcmpi(post.PReconMode, 'pointwise-exp')
            method = post.directThetaInterp;

            nFine = local_reconstruct_theta_matrix( ...
                thetaRaw, nRaw, thetaPlot, method);

            nFine = real(nFine);

            if post.clipNegative
                nFine(nFine < 0) = 0;
            end

            PPlot = trapz(x(:), nFine, 1).';
            PPlot = real(PPlot(:));
            PPlot(PPlot < 0) = 0;
            logPPlot = log(max(PPlot, realmin));
        else
            [PPlot, logPPlot, fitInfo] = local_reconstruct_from_logP( ...
                thetaRaw, logPRaw, thetaPlot, post, massForScaling, massSource);
        end

    case 'wkb'
        if ~isfield(S, 'W') || isempty(S.W) || ~isfield(S, 'u') || isempty(S.u)
            error('WKB 数据中需要同时包含 W 和 u。');
        end

        WRaw = local_to_xt_matrix(S.W, Nx, NthetaRaw, 'W');
        uRaw = local_extract_u_vector(S.u, Nx, NthetaRaw);

        epsVal = local_get_eps(S);

        Wpos = real(WRaw);
        Wpos(Wpos < 0) = 0;

        Araw = trapz(x(:), Wpos, 1).';
        Araw = real(Araw(:));
        Araw(Araw < realmin) = realmin;

        logPRaw = real(uRaw(:)) / epsVal + log(Araw);
        PRaw = local_exp_safe(logPRaw);
        [massForScaling, massSource] = local_mass_for_scaling(S, thetaRaw, logPRaw, post);

        if strcmpi(post.PReconMode, 'pointwise-exp')
            method = post.wkbThetaInterp;

            WFine = local_reconstruct_theta_matrix( ...
                thetaRaw, WRaw, thetaPlot, method);

            uFine = local_periodic_reconstruct( ...
                thetaRaw, uRaw, thetaPlot, method);

            WFine = real(WFine);

            if post.clipNegative
                WFine(WFine < 0) = 0;
            end

            expo = repmat(real(uFine(:).') / epsVal, Nx, 1);
            expo = min(max(expo, -700), 700);

            nFine = WFine .* exp(expo);

            PPlot = trapz(x(:), nFine, 1).';
            PPlot = real(PPlot(:));
            PPlot(PPlot < 0) = 0;
            logPPlot = log(max(PPlot, realmin));
        else
            [PPlot, logPPlot, fitInfo] = local_reconstruct_from_logP( ...
                thetaRaw, logPRaw, thetaPlot, post, massForScaling, massSource);
        end

    otherwise
        error('未知 solver 类型：%s。', solver);
end

massRaw = local_periodic_integral_log(thetaRaw, logPRaw);
massPlot = local_periodic_integral(thetaPlot, PPlot);

[one.maxPRaw, idRaw] = max(PRaw);
[one.maxPPlot, idPlot] = max(PPlot);

one.thetaRaw = thetaRaw(:);
one.thetaPlot = thetaPlot(:);

one.PRaw = PRaw(:);
one.PPlot = PPlot(:);

one.logPRaw = logPRaw(:);
one.logPPlot = logPPlot(:);

one.epsVal = epsVal;
one.Nx = Nx;
one.NthetaRaw = NthetaRaw;

one.massRaw = massRaw;
if exist('massForScaling', 'var')
    one.massForScaling = massForScaling;
    one.massSource = massSource;
end
one.massPlot = massPlot;

one.thetaPeakRaw = thetaRaw(idRaw);
one.thetaPeakPlot = thetaPlot(idPlot);

fields = fieldnames(fitInfo);
for i = 1:numel(fields)
    one.(fields{i}) = fitInfo.(fields{i});
end
end

function [PPlot, logPPlot, fitInfo] = local_reconstruct_from_logP(thetaRaw, logPRaw, thetaPlot, post, massForScaling, massSource)

fitInfo = struct();
if nargin < 5 || isempty(massForScaling)
    massForScaling = local_periodic_integral_log(thetaRaw, logPRaw);
end
if nargin < 6 || isempty(massSource)
    massSource = 'raw-data';
end
fitInfo.massForScaling = massForScaling;
fitInfo.massSource = massSource;

thetaRaw = thetaRaw(:);
logPRaw = real(logPRaw(:));
thetaPlot = thetaPlot(:);

switch lower(post.PReconMode)

    case 'log-pchip'
        logPPlot = local_periodic_interp1( ...
            thetaRaw, logPRaw, thetaPlot, post.logInterpMethod);

        massRaw = massForScaling;
        PPlot = local_shape_to_P(thetaPlot, logPPlot, massRaw, post.PMassMode);

    case 'log-quadratic-mass'
        [logPPlot, fitInfo] = local_log_quadratic_peak_fit( ...
            thetaRaw, logPRaw, thetaPlot, post);

        massRaw = massForScaling;
        PPlot = local_shape_to_P(thetaPlot, logPPlot, massRaw, post.PMassMode);

    case 'weno-anchored-quadratic-mass'
        [logPPlot, fitInfo] = local_weno_anchored_quadratic_peak_fit( ...
            thetaRaw, logPRaw, thetaPlot, post);

        massRaw = massForScaling;
        PPlot = local_shape_to_P(thetaPlot, logPPlot, massRaw, post.PMassMode);

    otherwise
        error('未知 PReconMode：%s。', post.PReconMode);
end

logPPlot = real(logPPlot(:));
PPlot = real(PPlot(:));
PPlot(PPlot < 0) = 0;
end

function [logPPlot, fitInfo] = local_weno_anchored_quadratic_peak_fit(thetaRaw, logPRaw, thetaPlot, post)

thetaRaw = thetaRaw(:);
logPRaw = real(logPRaw(:));
thetaPlot = thetaPlot(:);

K = numel(thetaRaw);

[thetaStar, ~, locInfo] = local_find_peak_location_from_logP( ...
    thetaRaw, logPRaw, post);


M = max(3, round(post.peakFitPoints));
M = min(M, K);

if mod(M, 2) == 0
    M = max(3, M - 1);
end

fitIndex = local_nearest_periodic_indices(thetaRaw, thetaStar, M);

zFit = local_periodic_distance(thetaRaw(fitIndex), thetaStar);
yFit = logPRaw(fitIndex);

[zFit, order] = sort(zFit);
yFit = yFit(order);
fitIndex = fitIndex(order);

yBase = max(yFit);
yShift = yFit - yBase;

if isfield(post, 'peakFitWeighted') && post.peakFitWeighted
    wFit = local_fit_weights(zFit, yShift, post);
else
    wFit = ones(size(zFit));
end

% 约束二次拟合：logP - yBase ≈ c2 z^2 + c0
X = [zFit.^2, ones(size(zFit))];
coef = local_weighted_linear_solve(X, yShift, wFit);

c2 = coef(1);
c0 = coef(2);

if c2 >= 0 || ~isfinite(c2)
    warning('anchored quadratic 拟合曲率非负，自动退回 log-quadratic-mass。');

    [logPPlot, fitInfo] = local_log_quadratic_peak_fit( ...
        thetaRaw, logPRaw, thetaPlot, post);

    fitInfo.fitFallback = 'log-quadratic-mass';
    return;
end

zPlot = local_periodic_distance(thetaPlot, thetaStar);
logPPlot = c2*zPlot.^2 + c0 + yBase;

fitInfo.fitIndex = fitIndex(:);
fitInfo.fitThetaCenter = thetaStar;
fitInfo.fitThetaPeak = thetaStar;
fitInfo.fitCurvature = c2;
fitInfo.fitInterceptShift = c0;
fitInfo.fitWeight = wFit(:);
fitInfo.fitFallback = 'none';
fitInfo.peakLocationMethod = post.peakLocationMethod;
fitInfo.peakLocationThetaRaw = locInfo.thetaRawPeak;
fitInfo.peakLocationThetaInterp = locInfo.thetaInterpPeak;
end

function [thetaStar, logPStar, info] = local_find_peak_location_from_logP(thetaRaw, logPRaw, post)

thetaRaw = thetaRaw(:);
logPRaw = real(logPRaw(:));

[~, idRaw] = max(logPRaw);
thetaRawPeak = thetaRaw(idRaw);

method = lower(post.peakLocationMethod);

if strcmpi(method, 'raw')
    thetaStar = thetaRawPeak;
    logPStar = logPRaw(idRaw);

    info.thetaRawPeak = thetaRawPeak;
    info.thetaInterpPeak = thetaStar;
    return;
end

Nsearch = max(post.peakSearchNtheta, 8*numel(thetaRaw));
thetaSearch = linspace(0, 1, Nsearch + 1).';
thetaSearch(end) = [];

logPSearch = local_periodic_reconstruct(thetaRaw, logPRaw, thetaSearch, method);

[logPStar, idMax] = max(logPSearch);
thetaStar = thetaSearch(idMax);

if isfield(post, 'peakRefineParabolic') && post.peakRefineParabolic
    thetaStar = local_parabolic_peak_refine(thetaSearch, logPSearch, idMax);
end

info.thetaRawPeak = thetaRawPeak;
info.thetaInterpPeak = thetaStar;
end

function thetaStar = local_parabolic_peak_refine(theta, y, idMax)

theta = theta(:);
y = real(y(:));

K = numel(theta);

if K < 3
    thetaStar = theta(idMax);
    return;
end

h = median(diff([theta; theta(1)+1]));

im = mod(idMax - 2, K) + 1;
i0 = idMax;
ip = mod(idMax, K) + 1;

ym = y(im);
y0 = y(i0);
yp = y(ip);

den = ym - 2*y0 + yp;

if abs(den) <= eps || ~isfinite(den)
    thetaStar = theta(i0);
    return;
end

delta = 0.5 * (ym - yp) / den;

if ~isfinite(delta) || abs(delta) > 1
    delta = 0;
end

thetaStar = mod(theta(i0) + delta*h, 1);
end

function [logPPlot, fitInfo] = local_log_quadratic_peak_fit(thetaRaw, logPRaw, thetaPlot, post)

thetaRaw = thetaRaw(:);
logPRaw = real(logPRaw(:));
thetaPlot = thetaPlot(:);

K = numel(thetaRaw);

[logPMaxRaw, kPeak] = max(logPRaw);
thetaCenter = thetaRaw(kPeak);

M = max(3, round(post.peakFitPoints));
M = min(M, K);

if mod(M, 2) == 0
    M = max(3, M - 1);
end

fitIndex = local_peak_window_indices(kPeak, K, M);

zFit = local_periodic_distance(thetaRaw(fitIndex), thetaCenter);
yFit = logPRaw(fitIndex);

[zFit, order] = sort(zFit);
yFit = yFit(order);
fitIndex = fitIndex(order);

yShift = yFit - logPMaxRaw;

deg = min(post.peakFitDegree, numel(zFit) - 1);
deg = max(1, deg);

if isfield(post, 'peakFitWeighted') && post.peakFitWeighted
    wFit = local_fit_weights(zFit, yShift, post);
    coefShift = local_weighted_polyfit(zFit, yShift, deg, wFit);
else
    coefShift = polyfit(zFit, yShift, deg);
    wFit = ones(size(zFit));
end

zPlot = local_periodic_distance(thetaPlot, thetaCenter);
logPPlot = polyval(coefShift, zPlot) + logPMaxRaw;

fitInfo.fitIndex = fitIndex(:);
fitInfo.fitThetaCenter = thetaCenter;
fitInfo.fitDegree = deg;
fitInfo.fitCoef = coefShift;
fitInfo.fitWeight = wFit(:);

if deg == 2
    fitInfo.fitCurvature = coefShift(1);

    if coefShift(1) < 0
        zStar = -coefShift(2) / (2 * coefShift(1));
        zMaxAllowed = max(abs(zFit));

        if abs(zStar) <= zMaxAllowed
            fitInfo.fitThetaPeak = mod(thetaCenter + zStar, 1);
        else
            fitInfo.fitThetaPeak = thetaCenter;
        end
    else
        fitInfo.fitThetaPeak = thetaCenter;
    end
else
    fitInfo.fitCurvature = NaN;
    fitInfo.fitThetaPeak = thetaCenter;
end

if deg == 2 && coefShift(1) >= 0
    warning('log-quadratic 拟合开口向上，自动退回 log-pchip。');

    logPPlot = local_periodic_interp1( ...
        thetaRaw, logPRaw, thetaPlot, post.logInterpMethod);

    fitInfo.fitFallback = 'log-pchip';
else
    fitInfo.fitFallback = 'none';
end
end

function PPlot = local_shape_to_P(thetaPlot, logShape, massRaw, massMode)

thetaPlot = thetaPlot(:);
logShape = real(logShape(:));

switch lower(massMode)

    case 'none'
        PPlot = local_exp_safe(logShape);

    case 'raw-data'
        logMax = max(logShape);
        shape = exp(logShape - logMax);
        intShape = local_periodic_integral(thetaPlot, shape);

        if intShape <= realmin || ~isfinite(intShape) || ~isfinite(massRaw)
            warning('raw-data 质量定标失败，退回 PMassMode = none。');
            PPlot = local_exp_safe(logShape);
        else
            PPlot = massRaw * shape / intShape;
        end

    case 'unit'
        logMax = max(logShape);
        shape = exp(logShape - logMax);
        intShape = local_periodic_integral(thetaPlot, shape);

        if intShape <= realmin || ~isfinite(intShape)
            PPlot = shape;
        else
            PPlot = shape / intShape;
        end

    otherwise
        error('未知 PMassMode：%s。', massMode);
end

PPlot = real(PPlot(:));
end

function coef = local_weighted_polyfit(x, y, deg, w)

x = x(:);
y = y(:);
w = w(:);

V = zeros(numel(x), deg + 1);

for j = 0:deg
    V(:, deg + 1 - j) = x.^j;
end

coef = local_weighted_linear_solve(V, y, w).';
end

function coef = local_weighted_linear_solve(X, y, w)

X = real(X);
y = real(y(:));
w = real(w(:));

w(~isfinite(w)) = 0;
w = max(w, 0);

if sum(w) <= realmin
    coef = X \ y;
    return;
end

sqrtw = sqrt(w);
XW = X .* sqrtw;
yW = y .* sqrtw;

coef = XW \ yW;
end

function w = local_fit_weights(zFit, yShift, post)

zFit = zFit(:);
yShift = yShift(:);

if ~isfield(post, 'peakFitWeightMode') || isempty(post.peakFitWeightMode)
    mode = 'log';
else
    mode = lower(post.peakFitWeightMode);
end

switch mode

    case 'none'
        w = ones(size(zFit));

    case 'log'
        if isfield(post, 'peakFitLogDrop') && ~isempty(post.peakFitLogDrop)
            drop = post.peakFitLogDrop;
        else
            drop = 8;
        end

        w = exp(max(yShift, -drop));
        w = max(w, exp(-drop));

    case 'gaussian'
        zScale = max(abs(zFit));

        if zScale <= realmin
            w = ones(size(zFit));
        else
            w = exp(-(zFit / (0.65*zScale)).^2);
        end

    otherwise
        error('未知 peakFitWeightMode：%s。', mode);
end

w = real(w(:));
end

function [x, theta] = local_get_grid(S)

if isfield(S, 'op') && isstruct(S.op) ...
        && isfield(S.op, 'x') && isfield(S.op, 'theta')
    x = S.op.x(:);
    theta = S.op.theta(:);
    return;
end

if isfield(S, 'grid') && isstruct(S.grid) ...
        && isfield(S.grid, 'x') && isfield(S.grid, 'theta')
    x = S.grid.x(:);
    theta = S.grid.theta(:);
    return;
end

if isfield(S, 'x') && isfield(S, 'theta')
    x = S.x(:);
    theta = S.theta(:);
    return;
end

if isfield(S, 'n') && ~isempty(S.n)
    sz = size(S.n);
elseif isfield(S, 'W') && ~isempty(S.W)
    sz = size(S.W);
else
    error('无法从文件中判断 x 和 theta 网格。');
end

x = linspace(0, 1, sz(1)).';
theta = linspace(0, 1, sz(2) + 1).';
theta(end) = [];
end

function A = local_to_xt_matrix(raw, Nx, Ntheta, name)

A = real(raw);

if isscalar(A)
    A = A * ones(Nx, Ntheta);

elseif isvector(A)
    v = A(:);

    if numel(v) == Nx
        A = repmat(v, 1, Ntheta);
    elseif numel(v) == Ntheta
        A = repmat(v(:).', Nx, 1);
    else
        error('%s 的向量长度不对，长度为 %d，但 Nx=%d, Ntheta=%d。', ...
            name, numel(v), Nx, Ntheta);
    end

elseif isequal(size(A), [Ntheta, Nx])
    A = A.';
end

if ~isequal(size(A), [Nx, Ntheta])
    error('%s 的矩阵大小为 %s，但期望大小是 [%d,%d]。', ...
        name, mat2str(size(A)), Nx, Ntheta);
end
end

function u = local_extract_u_vector(raw, Nx, Ntheta)

if isvector(raw)
    u = real(raw(:));
else
    U = local_to_xt_matrix(raw, Nx, Ntheta, 'u');
    u = real(U(1, :).');
end

if numel(u) ~= Ntheta
    error('u 的长度为 %d，但期望 Ntheta=%d。', numel(u), Ntheta);
end
end

function epsVal = local_get_eps(S)

if isfield(S, 'par') && isstruct(S.par) ...
        && isfield(S.par, 'eps') && ~isempty(S.par.eps)
    epsVal = S.par.eps;
elseif isfield(S, 'eps') && ~isempty(S.eps)
    epsVal = S.eps;
else
    error('WKB 文件中找不到 eps。');
end

epsVal = double(epsVal(1));
end

function epsVal = local_get_eps_if_possible(S)

epsVal = NaN;

if isfield(S, 'par') && isstruct(S.par) ...
        && isfield(S.par, 'eps') && ~isempty(S.par.eps)
    epsVal = double(S.par.eps(1));
elseif isfield(S, 'eps') && ~isempty(S.eps)
    epsVal = double(S.eps(1));
end
end

function AFine = local_reconstruct_theta_matrix(theta, A, thetaFine, method)

theta = theta(:);
thetaFine = thetaFine(:);

Nx = size(A, 1);
Kfine = numel(thetaFine);

AFine = zeros(Nx, Kfine);

for j = 1:Nx
    f = A(j, :).';
    AFine(j, :) = local_periodic_reconstruct(theta, f, thetaFine, method).';
end
end

function fq = local_periodic_reconstruct(theta, f, thetaq, method)

switch lower(method)
    case 'weno5'
        fq = local_weno5_periodic_interp(theta, f, thetaq);

    case {'linear','pchip','spline','makima'}
        fq = local_periodic_interp1(theta, f, thetaq, method);

    case 'none'
        fq = f(:);

    otherwise
        error('未知 theta 重构方法：%s', method);
end
end

function fq = local_periodic_interp1(theta, f, thetaq, method)

theta = theta(:);
f = real(f(:));
thetaq = mod(thetaq(:), 1);

[theta, idx] = sort(mod(theta, 1));
f = f(idx);

[theta, uniq] = unique(theta, 'stable');
f = f(uniq);

thetaExt = [theta(end) - 1; theta; theta(1) + 1];
fExt = [f(end); f; f(1)];

try
    fq = interp1(thetaExt, fExt, thetaq, method);
catch
    fq = interp1(thetaExt, fExt, thetaq, 'pchip');
end

fq = real(fq(:));
end

function fq = local_weno5_periodic_interp(theta, f, thetaq)

theta = theta(:);
f = real(f(:));
thetaq = mod(thetaq(:), 1);

K = numel(theta);

if K < 5
    fq = local_periodic_interp1(theta, f, thetaq, 'pchip');
    return;
end

dtheta = theta(2) - theta(1);
fq = zeros(size(thetaq));

small = 1.0e-12;

for m = 1:numel(thetaq)

    s = (thetaq(m) - theta(1)) / dtheta;
    i0 = floor(s) + 1;
    a = s - floor(s);

    im2 = local_pindex(i0 - 2, K);
    im1 = local_pindex(i0 - 1, K);
    i00 = local_pindex(i0,     K);
    ip1 = local_pindex(i0 + 1, K);
    ip2 = local_pindex(i0 + 2, K);

    fm2 = f(im2);
    fm1 = f(im1);
    f00 = f(i00);
    fp1 = f(ip1);
    fp2 = f(ip2);

    q0 = ((a + 1) * a / 2) * fm2 ...
       - ((a + 2) * a) * fm1 ...
       + ((a + 2) * (a + 1) / 2) * f00;

    q1 = (a * (a - 1) / 2) * fm1 ...
       - ((a + 1) * (a - 1)) * f00 ...
       + ((a + 1) * a / 2) * fp1;

    q2 = ((a - 1) * (a - 2) / 2) * f00 ...
       - (a * (a - 2)) * fp1 ...
       + (a * (a - 1) / 2) * fp2;

    d0 = a^2 / 12 - a / 4 + 1 / 6;
    d1 = 2 / 3 - a^2 / 6;
    d2 = a^2 / 12 + a / 4 + 1 / 6;

    d0 = max(d0, 0);
    d1 = max(d1, 0);
    d2 = max(d2, 0);

    dsum = d0 + d1 + d2;

    if dsum <= realmin
        fq(m) = local_periodic_interp1(theta, f, thetaq(m), 'pchip');
        continue;
    end

    d0 = d0 / dsum;
    d1 = d1 / dsum;
    d2 = d2 / dsum;

    beta0 = 13/12 * (fm2 - 2*fm1 + f00)^2 ...
          +  1/4 * (fm2 - 4*fm1 + 3*f00)^2;

    beta1 = 13/12 * (fm1 - 2*f00 + fp1)^2 ...
          +  1/4 * (fm1 - fp1)^2;

    beta2 = 13/12 * (f00 - 2*fp1 + fp2)^2 ...
          +  1/4 * (3*f00 - 4*fp1 + fp2)^2;

    alpha0 = d0 / (small + beta0)^2;
    alpha1 = d1 / (small + beta1)^2;
    alpha2 = d2 / (small + beta2)^2;

    asum = alpha0 + alpha1 + alpha2;

    if asum <= realmin || ~isfinite(asum)
        fq(m) = local_periodic_interp1(theta, f, thetaq(m), 'pchip');
    else
        w0 = alpha0 / asum;
        w1 = alpha1 / asum;
        w2 = alpha2 / asum;

        fq(m) = w0*q0 + w1*q1 + w2*q2;
    end
end

fq = real(fq(:));
end

function id = local_pindex(i, K)
id = mod(i - 1, K) + 1;
end

function idx = local_peak_window_indices(kPeak, K, M)

if M >= K
    idx = (1:K).';
    return;
end

r = floor(M / 2);
offsets = (-r:r).';
idx = mod(kPeak - 1 + offsets, K) + 1;
idx = idx(:);
end

function idx = local_nearest_periodic_indices(theta, theta0, M)

theta = theta(:);
d = abs(local_periodic_distance(theta, theta0));

[~, order] = sort(d, 'ascend');
idx = order(1:min(M, numel(theta)));
idx = sort(idx(:));
end

function d = local_periodic_distance(theta, theta0)
d = mod(theta(:) - theta0 + 0.5, 1) - 0.5;
end

function val = local_periodic_integral(theta, f)

theta = theta(:);
f = real(f(:));

K = numel(theta);

if K <= 1
    val = sum(f);
    return;
end

thetaSort = sort(mod(theta, 1));
dtheta = median(diff([thetaSort; thetaSort(1) + 1]));

val = dtheta * sum(f);
end

function val = local_periodic_integral_log(theta, logf)

theta = theta(:);
logf = real(logf(:));

K = numel(theta);

if K <= 1
    val = exp(min(max(logf(1), -700), 700));
    return;
end

thetaSort = sort(mod(theta, 1));
dtheta = median(diff([thetaSort; thetaSort(1) + 1]));

m = max(logf);

if ~isfinite(m)
    val = NaN;
    return;
end

if m > 700
    warning('logP 最大值超过 700，原始质量可能上溢。返回 Inf。');
    val = Inf;
else
    val = dtheta * exp(m) * sum(exp(logf - m));
end
end

function y = local_exp_safe(logy)
logy = real(logy(:));
y = exp(min(max(logy, -700), 700));
end

function [lineStyle, markerStyle] = local_plot_style(q, oneFile)

lineStyles = {'-', '--', '-.', ':', '-', '--', '-.', ':'};
markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x', '+'};

lineStyle = lineStyles{mod(q - 1, numel(lineStyles)) + 1};
markerStyle = markerStyles{mod(q - 1, numel(markerStyles)) + 1};

if isfield(oneFile, 'lineStyle') && ~isempty(oneFile.lineStyle)
    lineStyle = oneFile.lineStyle;
end

if isfield(oneFile, 'markerStyle') && ~isempty(oneFile.markerStyle)
    markerStyle = oneFile.markerStyle;
end
end

function markerIndex = local_select_raw_markers(thetaRaw, yRaw, maxRawMarkers, mode)

thetaRaw = thetaRaw(:);
yRaw = real(yRaw(:));

K = numel(thetaRaw);

if K == 0
    markerIndex = [];
    return;
end

if nargin < 3 || isempty(maxRawMarkers)
    maxRawMarkers = 8;
end

if nargin < 4 || isempty(mode)
    mode = 'uniform_peak';
end

maxRawMarkers = round(maxRawMarkers);

if maxRawMarkers <= 0
    markerIndex = [];
    return;
end

switch lower(mode)

    case 'all'
        markerIndex = (1:K).';
        return;

    case 'uniform'
        M = min(K, maxRawMarkers);
        markerIndex = unique(round(linspace(1, K, M))).';
        markerIndex = markerIndex(:);

    case 'uniform_peak'
        M = min(K, maxRawMarkers);

        if M == K
            markerIndex = (1:K).';
            return;
        end

        markerIndex = unique(round(linspace(1, K, M))).';
        markerIndex = markerIndex(:);

        [~, peakIndex] = max(yRaw);

        if ~ismember(peakIndex, markerIndex)
            if numel(markerIndex) >= M
                distToPeak = abs(markerIndex - peakIndex);
                [~, idReplace] = max(distToPeak);
                markerIndex(idReplace) = peakIndex;
            else
                markerIndex(end + 1) = peakIndex;
            end
        end

        markerIndex = unique(markerIndex(:));

        if numel(markerIndex) > M
            markerIndex = markerIndex(1:M);
        end

    otherwise
        error('未知 rawMarkerMode：%s。', mode);
end

markerIndex = markerIndex(markerIndex >= 1 & markerIndex <= K);
markerIndex = unique(markerIndex(:));
end