%% POST_COMPARE_WKB_RAW_U_W_N
% -------------------------------------------------------------------------
% 比较若干个 WKB 结果中的原始 u, W, n。
%
% 注意：
%   1. 本脚本不做 theta 方向插值。
%   2. Ntheta 不同的数据，各自画在自己的原始 theta 网格上。
%   3. W 和 n 是 x-theta 二维量，这里默认取一个 x 切片比较。
%   4. 同时额外画 P(theta)=int_x n(x,theta) dx，便于检查 trait 分布。
% -------------------------------------------------------------------------

clear; clc; close all;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

caseName = 'test3_wkb_advantage_refinement';
caseDir = fullfile(root, 'data', caseName);

% -------------------------------------------------------------------------
% 手动指定要比较的数据。
% -------------------------------------------------------------------------

files = { ...
    struct( ...
        'dirName',    'wkb_eps0p0001_dt0p0001_Nx40_Ntheta15_t2', ...
        'fileName',   'auto', ...
        'targetTime', 2, ...
        'solver',     'wkb', ...
        'label',      'WKB, N_\theta=15'), ...
    struct( ...
        'dirName',    'wkb_eps0p0001_dt0p0001_Nx40_Ntheta16_t2', ...
        'fileName',   'auto', ...
        'targetTime', 2, ...
        'solver',     'wkb', ...
        'label',      'WKB, N_\theta=16'), ...
    struct( ...
        'dirName',    'wkb_eps0p0001_dt0p0001_Nx40_Ntheta17_t2', ...
        'fileName',   'auto', ...
        'targetTime', 2, ...
        'solver',     'wkb', ...
        'label',      'WKB, N_\theta=17') ...
    };

% -------------------------------------------------------------------------
% 比较设置
% -------------------------------------------------------------------------

post.xSliceMode = 'middle';      % {'middle','index','nearest_x'}
post.xIndex = 20;                % xSliceMode='index' 时使用
post.xValue = 0.5;               % xSliceMode='nearest_x' 时使用

post.uDisplay = 'u_over_eps';    % {'u','u_over_eps'}
post.nDisplay = 'raw';    % {'raw','normalized','log'}
post.PDisplay = 'raw';    % {'raw','normalized','log'}

post.showMarkers = true;
post.markerSize = 6.5;
post.lineWidth = 2.0;

post.saveFigures = true;
post.outputDir = fullfile(root, 'data', 'figures', caseName);
post.figurePrefix = 'compare_wkb_raw_u_W_n';

if exist(post.outputDir, 'dir') ~= 7
    mkdir(post.outputDir);
end

%% ========================= 读取数据 =========================

data = cell(numel(files), 1);

for q = 1:numel(files)

    dataDir = fullfile(caseDir, files{q}.dirName);

    if strcmpi(files{q}.fileName, 'auto')
        filePath = local_find_snapshot_by_time( ...
            dataDir, files{q}.solver, files{q}.targetTime);
    else
        filePath = fullfile(dataDir, files{q}.fileName);
        if exist(filePath, 'file') ~= 2
            error('找不到文件：\n%s', filePath);
        end
    end

    Sraw = load(filePath);
    S = local_select_data_struct(Sraw);

    one = local_extract_wkb_raw(S, post);
    one.label = files{q}.label;
    one.filePath = filePath;

    data{q} = one;

    fprintf('\n[%s]\n', one.label);
    fprintf('  file     = %s\n', filePath);
    fprintf('  eps      = %.16g\n', one.epsVal);
    fprintf('  Nx       = %d\n', one.Nx);
    fprintf('  Ntheta   = %d\n', one.Ntheta);
    fprintf('  xIndex   = %d\n', one.xIndex);
    fprintf('  xValue   = %.16g\n', one.xValue);
    fprintf('  u peak   = theta %.8f\n', one.thetaUPeak);
    fprintf('  n peak   = theta %.8f at selected x\n', one.thetaNPeak);
    fprintf('  P peak   = theta %.8f\n', one.thetaPPeak);
end

%% ========================= 图 1：u, W, n, P 原始曲线 =========================

fig1 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [3, 3, 24, 16]);

tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% u(theta)
ax1 = nexttile;
hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    switch lower(post.uDisplay)
        case 'u'
            y = one.u(:);
            ylab = 'u(\theta)';
        case 'u_over_eps'
            y = one.u(:) / one.epsVal;
            ylab = 'u(\theta)/\epsilon';
        otherwise
            error('未知 uDisplay：%s', post.uDisplay);
    end

    local_plot_raw_curve(ax1, one.theta, y, one.label, ls, mk, post);
end

xlabel(ax1, '\theta');
ylabel(ax1, ylab);
title(ax1, 'phase variable');

% W(x*,theta)
ax2 = nexttile;
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    y = one.Wslice(:);
    local_plot_raw_curve(ax2, one.theta, y, one.label, ls, mk, post);
end

xlabel(ax2, '\theta');
ylabel(ax2, 'W(x_\ast,\theta)');
title(ax2, sprintf('amplitude at selected x'));

% n(x*,theta)
ax3 = nexttile;
hold(ax3, 'on'); box(ax3, 'on'); grid(ax3, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    y = local_display_transform(one.nSlice(:), post.nDisplay);
    local_plot_raw_curve(ax3, one.theta, y, one.label, ls, mk, post);
end

xlabel(ax3, '\theta');

switch lower(post.nDisplay)
    case 'raw'
        ylabel(ax3, 'n(x_\ast,\theta)');
    case 'normalized'
        ylabel(ax3, 'normalized n(x_\ast,\theta)');
    case 'log'
        ylabel(ax3, 'log_{10}(n(x_\ast,\theta))');
end

title(ax3, 'reconstructed density at selected x');

% P(theta)
ax4 = nexttile;
hold(ax4, 'on'); box(ax4, 'on'); grid(ax4, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    y = local_display_transform(one.Ptheta(:), post.PDisplay);
    local_plot_raw_curve(ax4, one.theta, y, one.label, ls, mk, post);
end

xlabel(ax4, '\theta');

switch lower(post.PDisplay)
    case 'raw'
        ylabel(ax4, 'P(\theta)');
    case 'normalized'
        ylabel(ax4, 'normalized P(\theta)');
    case 'log'
        ylabel(ax4, 'log_{10}(P(\theta))');
end

title(ax4, 'trait marginal');

legend(ax4, 'Location', 'best', 'Box', 'off');

for ax = [ax1, ax2, ax3, ax4]
    xlim(ax, [0, 1]);
    set(ax, 'FontSize', 11);
end

if post.saveFigures
    pngFile = fullfile(post.outputDir, [post.figurePrefix, '_curves.png']);
    figFile = fullfile(post.outputDir, [post.figurePrefix, '_curves.fig']);

    try
        exportgraphics(fig1, pngFile, 'Resolution', 300);
    catch
        saveas(fig1, pngFile);
    end

    savefig(fig1, figFile);

    fprintf('\n曲线图已保存：\n%s\n%s\n', pngFile, figFile);
end

%% ========================= 图 2：W 和 n 的二维热图 =========================

fig2 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [4, 4, 24, 6*numel(data)]);

tiledlayout(fig2, numel(data), 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for q = 1:numel(data)

    one = data{q};

    axW = nexttile;
    imagesc(axW, one.theta, one.x, one.W);
    set(axW, 'YDir', 'normal');
    colorbar(axW);
    xlabel(axW, '\theta');
    ylabel(axW, 'x');
    title(axW, ['W, ', one.label], 'Interpreter', 'tex');
    set(axW, 'FontSize', 11);

    axN = nexttile;

    switch lower(post.nDisplay)
        case 'raw'
            nShow = one.n;
            nTitle = 'n';
        case 'normalized'
            nMax = max(one.n(:));
            if isfinite(nMax) && nMax > 0
                nShow = one.n / nMax;
            else
                nShow = one.n;
            end
            nTitle = 'normalized n';
        case 'log'
            nShow = log10(max(one.n, realmin));
            nTitle = 'log_{10} n';
        otherwise
            error('未知 nDisplay：%s', post.nDisplay);
    end

    imagesc(axN, one.theta, one.x, nShow);
    set(axN, 'YDir', 'normal');
    colorbar(axN);
    xlabel(axN, '\theta');
    ylabel(axN, 'x');
    title(axN, [nTitle, ', ', one.label], 'Interpreter', 'tex');
    set(axN, 'FontSize', 11);
end

if post.saveFigures
    pngFile = fullfile(post.outputDir, [post.figurePrefix, '_heatmaps.png']);
    figFile = fullfile(post.outputDir, [post.figurePrefix, '_heatmaps.fig']);

    try
        exportgraphics(fig2, pngFile, 'Resolution', 300);
    catch
        saveas(fig2, pngFile);
    end

    savefig(fig2, figFile);

    fprintf('\n热图已保存：\n%s\n%s\n', pngFile, figFile);
end

%% ========================================================================
% 局部函数
% ========================================================================

function filePath = local_find_snapshot_by_time(dataDir, solver, targetTime)
% 在 dataDir 中查找最接近 targetTime 的 snapshot 文件。

if exist(dataDir, 'dir') ~= 7
    error('数据目录不存在：\n%s', dataDir);
end

pat = sprintf('snapshot_%s*.mat', lower(solver));
d = dir(fullfile(dataDir, pat));

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
end

function t = local_parse_time_from_snapshot_name(name)
% 从文件名中解析 t。
%
% 支持：
%   snapshot_wkb_t2.mat
%   snapshot_wkb_t0p5.mat
%   snapshot_wkb_t2p5em4.mat
%   snapshot_wkb_0002000_t2.mat
%   snapshot_wkb_step0002000_t2.mat

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
% 兼容不同保存格式。

S = L;

if isfield(L, 'snapshot') && isstruct(L.snapshot)
    S = L.snapshot;
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

function one = local_extract_wkb_raw(S, post)
% 从一个 WKB snapshot 中取出原始 u, W，并重构 n。

[x, theta] = local_get_grid(S);

Nx = numel(x);
Ntheta = numel(theta);

if ~isfield(S, 'W') || isempty(S.W)
    error('WKB 数据中没有 W。');
end

if ~isfield(S, 'u') || isempty(S.u)
    error('WKB 数据中没有 u。');
end

W = local_to_xt_matrix(S.W, Nx, Ntheta, 'W');

uRaw = S.u;
if isvector(uRaw)
    u = real(uRaw(:));
    U = repmat(u(:).', Nx, 1);
else
    U = local_to_xt_matrix(uRaw, Nx, Ntheta, 'u');
    u = real(U(1, :).');
end

epsVal = local_get_eps(S);

expo = real(U) / epsVal;
expo = min(max(expo, -700), 700);

Wpos = real(W);
Wpos(Wpos < 0) = 0;

n = Wpos .* exp(expo);

xIndex = local_choose_x_index(x, post);

Wslice = W(xIndex, :).';
nSlice = n(xIndex, :).';

Ptheta = trapz(x(:), n, 1).';
Ptheta = real(Ptheta(:));

[~, idU] = max(u);
[~, idN] = max(nSlice);
[~, idP] = max(Ptheta);

one.x = x(:);
one.theta = theta(:);
one.W = real(W);
one.u = real(u(:));
one.n = real(n);
one.Wslice = real(Wslice(:));
one.nSlice = real(nSlice(:));
one.Ptheta = real(Ptheta(:));

one.epsVal = epsVal;
one.Nx = Nx;
one.Ntheta = Ntheta;
one.xIndex = xIndex;
one.xValue = x(xIndex);

one.thetaUPeak = theta(idU);
one.thetaNPeak = theta(idN);
one.thetaPPeak = theta(idP);
end

function xIndex = local_choose_x_index(x, post)
% 选择 W 和 n 的 x 切片位置。

Nx = numel(x);

switch lower(post.xSliceMode)

    case 'middle'
        xIndex = round((Nx + 1) / 2);

    case 'index'
        xIndex = post.xIndex;
        xIndex = max(1, min(Nx, round(xIndex)));

    case 'nearest_x'
        [~, xIndex] = min(abs(x(:) - post.xValue));

    otherwise
        error('未知 xSliceMode：%s。可选：middle, index, nearest_x。', post.xSliceMode);
end
end

function [x, theta] = local_get_grid(S)
% 从数据结构中读取 x 和 theta。

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

if isfield(S, 'W') && ~isempty(S.W)
    sz = size(S.W);
else
    error('无法从文件中判断 x 和 theta 网格。');
end

x = linspace(0, 1, sz(1)).';
theta = linspace(0, 1, sz(2) + 1).';
theta(end) = [];
end

function A = local_to_xt_matrix(raw, Nx, Ntheta, name)
% 统一转成大小为 Nx × Ntheta 的矩阵。

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

function epsVal = local_get_eps(S)
% 从数据结构中读取 epsilon。

if isfield(S, 'par') && isstruct(S.par) ...
        && isfield(S.par, 'eps') && ~isempty(S.par.eps)
    epsVal = S.par.eps;
elseif isfield(S, 'eps') && ~isempty(S.eps)
    epsVal = S.eps;
else
    error('WKB 文件中找不到 eps。请检查保存数据是否包含 par.eps 或 eps。');
end

epsVal = double(epsVal(1));
end

function y = local_display_transform(yRaw, mode)
% 对 n 或 P 做显示变换，不改变原始数据。

yRaw = real(yRaw(:));

switch lower(mode)

    case 'raw'
        y = yRaw;

    case 'normalized'
        ymax = max(yRaw);
        if isfinite(ymax) && ymax > 0
            y = yRaw / ymax;
        else
            y = yRaw;
        end

    case 'log'
        y = log10(max(yRaw, realmin));

    otherwise
        error('未知显示模式：%s。可选：raw, normalized, log。', mode);
end
end

function local_plot_raw_curve(ax, theta, y, label, lineStyle, markerStyle, post)
% 画原始 theta 网格上的曲线，不做插值。

theta = theta(:);
y = y(:);

if post.showMarkers
    plot(ax, theta, y, ...
        'LineStyle', lineStyle, ...
        'Marker', markerStyle, ...
        'MarkerSize', post.markerSize, ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', label);
else
    plot(ax, theta, y, ...
        'LineStyle', lineStyle, ...
        'Marker', 'none', ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', label);
end
end

function [lineStyle, markerStyle] = local_plot_style(q)
% 不同曲线使用不同线型和 marker。

lineStyles = {'-', '--', '-.', ':', '-', '--', '-.', ':'};
markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x', '+'};

lineStyle = lineStyles{mod(q - 1, numel(lineStyles)) + 1};
markerStyle = markerStyles{mod(q - 1, numel(markerStyles)) + 1};
end