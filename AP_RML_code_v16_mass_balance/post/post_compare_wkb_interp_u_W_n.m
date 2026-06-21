%% POST_COMPARE_WKB_INTERP_U_W_N
% -------------------------------------------------------------------------
% 比较若干个 WKB 结果的插值后 u, W, n, P。
%
% 功能：
%   1. 对 u(theta) 做 theta 方向插值；
%   2. 对 W(x,theta) 做 theta 方向插值；
%   3. 用插值后的 W 和 u 重构
%          n(x,theta) = W(x,theta) * exp(u(theta)/eps);
%   4. 画插值后的
%          u(theta), W(x*,theta), n(x*,theta), P(theta)
%      同时叠加原始离散点，便于比较“原始数据”和“插值结果”。
%
% 注意：
%   1. 这里只处理 WKB 数据。
%   2. WKB 插值顺序固定为：先插值 u 和 W，再重构 n。
%   3. n 和 P 默认不做 normalize。
% -------------------------------------------------------------------------

clear; clc; close all;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

caseName = 'test3_wkb_advantage_refinement';
caseDir = fullfile(root, 'data', caseName);

% -------------------------------------------------------------------------
% 要比较的数据
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
% 插值与显示选项
% -------------------------------------------------------------------------
post.thetaInterpMethod = 'weno5';    % {'none','weno5','pchip','spline','linear','makima'}
post.NthetaPlot = 1024;              % 插值后细 theta 网格点数

post.xSliceMode = 'middle';          % {'middle','index','nearest_x'}
post.xIndex = 20;                    % xSliceMode='index' 时使用
post.xValue = 0.5;                   % xSliceMode='nearest_x' 时使用

post.uDisplay = 'u_over_eps';        % {'u','u_over_eps'}
post.nDisplay = 'raw';               % {'raw','normalized','log'}
post.PDisplay = 'raw';               % {'raw','normalized','log'}

post.clipNegativeW = true;           % 插值后 W 是否截断为非负
post.showRawPoints = true;           % 是否叠加原始离散点
post.maxRawMarkers = 8;              % 每条曲线最多显示几个原始点
post.rawMarkerMode = 'uniform_peak'; % {'uniform','uniform_peak','all'}

post.lineWidth = 2.3;
post.rawMarkerSize = 6.5;
post.rawMarkerLineWidth = 1.3;

post.saveFigures = true;
post.outputDir = fullfile(root, 'data', 'figures', caseName);
post.figurePrefix = 'compare_wkb_interp_u_W_n';

if exist(post.outputDir, 'dir') ~= 7
    mkdir(post.outputDir);
end

%% ========================= 读取并插值 =========================

data = cell(numel(files), 1);

fprintf('\n================ 插值后 WKB 比较 ================\n');
fprintf('theta 插值方式：%s\n', post.thetaInterpMethod);
fprintf('插值细网格点数：%d\n', post.NthetaPlot);

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

    one = local_extract_wkb_interp(S, post);
    one.label = files{q}.label;
    one.filePath = filePath;

    data{q} = one;

    fprintf('\n[%s]\n', one.label);
    fprintf('  file       = %s\n', filePath);
    fprintf('  eps        = %.16g\n', one.epsVal);
    fprintf('  Nx         = %d\n', one.Nx);
    fprintf('  NthetaRaw  = %d\n', one.NthetaRaw);
    fprintf('  xIndex     = %d\n', one.xIndex);
    fprintf('  xValue     = %.16g\n', one.xValue);
    fprintf('  raw u peak = %.8f\n', one.thetaUPeakRaw);
    fprintf('  int u peak = %.8f\n', one.thetaUPeakInterp);
    fprintf('  raw P peak = %.8f\n', one.thetaPPeakRaw);
    fprintf('  int P peak = %.8f\n', one.thetaPPeakInterp);
end

%% ========================= 图 1：插值后的 u, W, n, P =========================

fig1 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [3, 3, 24, 16]);

tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% -------------------------------------------------------------------------
% 1) u(theta)
% -------------------------------------------------------------------------
ax1 = nexttile;
hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    switch lower(post.uDisplay)
        case 'u'
            yPlot = one.uInterp(:);
            yRaw  = one.uRaw(:);
            ylab = 'u(\theta)';
        case 'u_over_eps'
            yPlot = one.uInterp(:) / one.epsVal;
            yRaw  = one.uRaw(:) / one.epsVal;
            ylab = 'u(\theta)/\epsilon';
        otherwise
            error('未知 uDisplay：%s', post.uDisplay);
    end

    hLine = plot(ax1, one.thetaPlot, yPlot, ...
        'LineStyle', ls, ...
        'Marker', 'none', ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', one.label);

    if post.showRawPoints
        idx = local_select_raw_markers(one.thetaRaw, yRaw, ...
            post.maxRawMarkers, post.rawMarkerMode);

        plot(ax1, one.thetaRaw(idx), yRaw(idx), ...
            'LineStyle', 'none', ...
            'Marker', mk, ...
            'MarkerSize', post.rawMarkerSize, ...
            'LineWidth', post.rawMarkerLineWidth, ...
            'Color', hLine.Color, ...
            'MarkerFaceColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

xlabel(ax1, '\theta');
ylabel(ax1, ylab);
title(ax1, 'interpolated phase variable');

% -------------------------------------------------------------------------
% 2) W(x*,theta)
% -------------------------------------------------------------------------
ax2 = nexttile;
hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    yPlot = one.WsliceInterp(:);
    yRaw  = one.WsliceRaw(:);

    hLine = plot(ax2, one.thetaPlot, yPlot, ...
        'LineStyle', ls, ...
        'Marker', 'none', ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', one.label);

    if post.showRawPoints
        idx = local_select_raw_markers(one.thetaRaw, yRaw, ...
            post.maxRawMarkers, post.rawMarkerMode);

        plot(ax2, one.thetaRaw(idx), yRaw(idx), ...
            'LineStyle', 'none', ...
            'Marker', mk, ...
            'MarkerSize', post.rawMarkerSize, ...
            'LineWidth', post.rawMarkerLineWidth, ...
            'Color', hLine.Color, ...
            'MarkerFaceColor', 'none', ...
            'HandleVisibility', 'off');
    end
end

xlabel(ax2, '\theta');
ylabel(ax2, 'W(x_\ast,\theta)');
title(ax2, 'interpolated amplitude at selected x');

% -------------------------------------------------------------------------
% 3) n(x*,theta)
% -------------------------------------------------------------------------
ax3 = nexttile;
hold(ax3, 'on'); box(ax3, 'on'); grid(ax3, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    yPlot = local_display_transform(one.nSliceInterp(:), post.nDisplay);
    yRaw  = local_display_transform(one.nSliceRaw(:), post.nDisplay);

    hLine = plot(ax3, one.thetaPlot, yPlot, ...
        'LineStyle', ls, ...
        'Marker', 'none', ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', one.label);

    if post.showRawPoints
        idx = local_select_raw_markers(one.thetaRaw, yRaw, ...
            post.maxRawMarkers, post.rawMarkerMode);

        plot(ax3, one.thetaRaw(idx), yRaw(idx), ...
            'LineStyle', 'none', ...
            'Marker', mk, ...
            'MarkerSize', post.rawMarkerSize, ...
            'LineWidth', post.rawMarkerLineWidth, ...
            'Color', hLine.Color, ...
            'MarkerFaceColor', 'none', ...
            'HandleVisibility', 'off');
    end
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

% -------------------------------------------------------------------------
% 4) P(theta)
% -------------------------------------------------------------------------
ax4 = nexttile;
hold(ax4, 'on'); box(ax4, 'on'); grid(ax4, 'on');

for q = 1:numel(data)
    one = data{q};
    [ls, mk] = local_plot_style(q);

    yPlot = local_display_transform(one.PthetaInterp(:), post.PDisplay);
    yRaw  = local_display_transform(one.PthetaRaw(:), post.PDisplay);

    hLine = plot(ax4, one.thetaPlot, yPlot, ...
        'LineStyle', ls, ...
        'Marker', 'none', ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', one.label);

    if post.showRawPoints
        idx = local_select_raw_markers(one.thetaRaw, yRaw, ...
            post.maxRawMarkers, post.rawMarkerMode);

        plot(ax4, one.thetaRaw(idx), yRaw(idx), ...
            'LineStyle', 'none', ...
            'Marker', mk, ...
            'MarkerSize', post.rawMarkerSize, ...
            'LineWidth', post.rawMarkerLineWidth, ...
            'Color', hLine.Color, ...
            'MarkerFaceColor', 'none', ...
            'HandleVisibility', 'off');
    end
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

title(ax4, 'interpolated trait marginal');

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

    fprintf('\n插值曲线图已保存：\n%s\n%s\n', pngFile, figFile);
end

%% ========================= 图 2：插值后的 u 和 W 热图辅助检查 =========================

fig2 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [4, 4, 24, 6*numel(data)]);

tiledlayout(fig2, numel(data), 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for q = 1:numel(data)
    one = data{q};

    % u(theta) 复制成 x-theta 图，仅辅助观察
    axU = nexttile;
    Ushow = repmat(one.uInterp(:).', one.Nx, 1);
    imagesc(axU, one.thetaPlot, one.x, Ushow);
    set(axU, 'YDir', 'normal');
    colorbar(axU);
    xlabel(axU, '\theta');
    ylabel(axU, 'x');
    title(axU, ['interpolated u, ', one.label], 'Interpreter', 'tex');
    set(axU, 'FontSize', 11);

    axW = nexttile;
    imagesc(axW, one.thetaPlot, one.x, one.WInterp);
    set(axW, 'YDir', 'normal');
    colorbar(axW);
    xlabel(axW, '\theta');
    ylabel(axW, 'x');
    title(axW, ['interpolated W, ', one.label], 'Interpreter', 'tex');
    set(axW, 'FontSize', 11);
end

if post.saveFigures
    pngFile = fullfile(post.outputDir, [post.figurePrefix, '_uW_heatmaps.png']);
    figFile = fullfile(post.outputDir, [post.figurePrefix, '_uW_heatmaps.fig']);

    try
        exportgraphics(fig2, pngFile, 'Resolution', 300);
    catch
        saveas(fig2, pngFile);
    end
    savefig(fig2, figFile);

    fprintf('\n插值后 u/W 热图已保存：\n%s\n%s\n', pngFile, figFile);
end

%% ========================================================================
% 局部函数
% ========================================================================

function filePath = local_find_snapshot_by_time(dataDir, solver, targetTime)
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

function one = local_extract_wkb_interp(S, post)
% 读取原始 WKB 数据，并在 theta 方向做插值，再重构 n。

[x, thetaRaw] = local_get_grid(S);

Nx = numel(x);
NthetaRaw = numel(thetaRaw);

if ~isfield(S, 'W') || isempty(S.W)
    error('WKB 数据中没有 W。');
end

if ~isfield(S, 'u') || isempty(S.u)
    error('WKB 数据中没有 u。');
end

WRaw = local_to_xt_matrix(S.W, Nx, NthetaRaw, 'W');

uRawData = S.u;
if isvector(uRawData)
    uRaw = real(uRawData(:));
    URawMat = repmat(uRaw(:).', Nx, 1);
else
    URawMat = local_to_xt_matrix(uRawData, Nx, NthetaRaw, 'u');
    uRaw = real(URawMat(1, :).');
end

epsVal = local_get_eps(S);

% 原始重构 n
expoRaw = URawMat / epsVal;
expoRaw = min(max(real(expoRaw), -700), 700);

WRawPos = real(WRaw);
WRawPos(WRawPos < 0) = 0;

nRaw = WRawPos .* exp(expoRaw);
PthetaRaw = trapz(x(:), nRaw, 1).';
PthetaRaw = real(PthetaRaw(:));

% 插值后网格
thetaPlot = linspace(0, 1, post.NthetaPlot + 1).';
thetaPlot(end) = [];

% 插值 u 和 W
if strcmpi(post.thetaInterpMethod, 'none')
    thetaPlot = thetaRaw;
    uInterp = uRaw;
    WInterp = WRaw;
else
    uInterp = local_periodic_reconstruct(thetaRaw, uRaw, thetaPlot, post.thetaInterpMethod);

    WInterp = local_reconstruct_theta_matrix( ...
        thetaRaw, WRaw, thetaPlot, post.thetaInterpMethod);
end

if post.clipNegativeW
    WInterp = real(WInterp);
    WInterp(WInterp < 0) = 0;
end

% 插值后重构 n
UInterpMat = repmat(real(uInterp(:).'), Nx, 1);
expoInterp = UInterpMat / epsVal;
expoInterp = min(max(real(expoInterp), -700), 700);

nInterp = WInterp .* exp(expoInterp);
PthetaInterp = trapz(x(:), nInterp, 1).';
PthetaInterp = real(PthetaInterp(:));

% 取 x 切片
xIndex = local_choose_x_index(x, post);

WsliceRaw = WRaw(xIndex, :).';
WsliceInterp = WInterp(xIndex, :).';

nSliceRaw = nRaw(xIndex, :).';
nSliceInterp = nInterp(xIndex, :).';

[~, idURaw] = max(uRaw);
[~, idUInt] = max(uInterp);
[~, idPRaw] = max(PthetaRaw);
[~, idPInt] = max(PthetaInterp);

one.x = x(:);
one.thetaRaw = thetaRaw(:);
one.thetaPlot = thetaPlot(:);

one.WRaw = real(WRaw);
one.WInterp = real(WInterp);

one.uRaw = real(uRaw(:));
one.uInterp = real(uInterp(:));

one.nRaw = real(nRaw);
one.nInterp = real(nInterp);

one.WsliceRaw = real(WsliceRaw(:));
one.WsliceInterp = real(WsliceInterp(:));

one.nSliceRaw = real(nSliceRaw(:));
one.nSliceInterp = real(nSliceInterp(:));

one.PthetaRaw = real(PthetaRaw(:));
one.PthetaInterp = real(PthetaInterp(:));

one.epsVal = epsVal;
one.Nx = Nx;
one.NthetaRaw = NthetaRaw;
one.xIndex = xIndex;
one.xValue = x(xIndex);

one.thetaUPeakRaw = thetaRaw(idURaw);
one.thetaUPeakInterp = thetaPlot(idUInt);
one.thetaPPeakRaw = thetaRaw(idPRaw);
one.thetaPPeakInterp = thetaPlot(idPInt);
end

function xIndex = local_choose_x_index(x, post)
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
        error('未知 theta 插值方法：%s', method);
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

thetaExt = [theta(end)-1; theta; theta(1)+1];
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

function y = local_display_transform(yRaw, mode)
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

function [lineStyle, markerStyle] = local_plot_style(q)
lineStyles = {'-', '--', '-.', ':', '-', '--', '-.', ':'};
markerStyles = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h', 'x', '+'};

lineStyle = lineStyles{mod(q - 1, numel(lineStyles)) + 1};
markerStyle = markerStyles{mod(q - 1, numel(markerStyles)) + 1};
end

function markerIndex = local_select_raw_markers(thetaRaw, yRaw, maxRawMarkers, mode)
thetaRaw = thetaRaw(:);
yRaw = yRaw(:);
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
        error('未知 rawMarkerMode：%s。可选：uniform, uniform_peak, all。', mode);
end

markerIndex = markerIndex(markerIndex >= 1 & markerIndex <= K);
markerIndex = unique(markerIndex(:));
end