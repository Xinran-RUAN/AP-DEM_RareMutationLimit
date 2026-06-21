%% POST_TEST3_PTHETA_SELECTED_FILES
% -------------------------------------------------------------------------
% 例 3 专用后处理：
%   手动指定若干个 direct / WKB 数据目录，画 P(theta) 对比图。
%
% direct:
%   对 n(x,theta) 在 theta 方向插值到细网格，
%   再计算 P(theta) = int_x n(x,theta) dx。
%
% WKB:
%   分别对 W(x,theta) 和 u(theta) 在 theta 方向插值到细网格，
%   再计算 n = W exp(u/eps)，
%   最后计算 P(theta) = int_x n(x,theta) dx。
%
% 注意：
%   1. 所有重构只在后处理中进行，不改核心求解器。
%   2. 插值后的光滑曲线画在 thetaPlot 上。
%   3. marker 画在插值前的原始位置 thetaRaw, PRaw 上。
% -------------------------------------------------------------------------

clear; clc; close all;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

caseName = 'test3_wkb_advantage_refinement';
caseDir = fullfile(root, 'data', caseName);

% -------------------------------------------------------------------------
% 后处理重构选项
% -------------------------------------------------------------------------
% directThetaInterp 控制 direct 结果中 n(x,theta) 的 theta 方向插值。
% wkbThetaInterp    控制 WKB 结果中 W(x,theta) 和 u(theta) 的 theta 方向插值。
%
% 可选：
%   'none'    不插值，直接画原始 theta 网格上的 P(theta)
%   'weno5'   WENO5 型后处理插值
%   'pchip'   pchip 周期插值
%   'spline'  spline 周期插值
%   'linear'  linear 周期插值
%   'makima'  makima 周期插值
% -------------------------------------------------------------------------

post.directThetaInterp = 'pchip';     % {'none','weno5','pchip','spline','linear','makima'}
post.wkbThetaInterp    = 'weno5';     % {'none','weno5','pchip','spline','linear','makima'}

% -------------------------------------------------------------------------
% 画图选项
% -------------------------------------------------------------------------
post.NthetaPlot = 4096;               % 后处理细 theta 网格点数，越大曲线越光滑

post.showRawPoints = true;            % 是否叠加插值前的原始离散点
post.maxRawMarkers = 8;               % 每条曲线最多画几个原始点 marker
post.rawMarkerMode = 'uniform_peak';  % {'uniform','uniform_peak','all'}

post.lineWidth = 2.3;                 % 光滑曲线线宽
post.rawMarkerSize = 6.5;             % 原始点 marker 大小
post.rawMarkerLineWidth = 1.3;        % 原始点 marker 线宽

post.normalizeP = false;               % 是否归一化 P(theta)
post.clipNegative = false;             % 重构后是否截断负值

post.outputDir = fullfile(root, 'data', 'figures', caseName);
post.figureName = 'test3_Ptheta_selected_files';

if exist(post.outputDir, 'dir') ~= 7
    mkdir(post.outputDir);
end

% -------------------------------------------------------------------------
% 手动指定需要比较的数据。
%
% dirName:
%   data/caseName 下的子文件夹名。
%
% fileName:
%   可以填具体文件名，例如：
%       'snapshot_direct_t2.mat'
%       'snapshot_wkb_t2.mat'
%       'result_direct_final.mat'
%       'result_wkb_final.mat'
%
%   也可以填 'auto'，此时程序会在对应目录中自动找最接近 targetTime
%   的 snapshot 文件。
%
% solver:
%   'direct' 或 'wkb'。
%
% label:
%   图例名称。
%
% reconMethod:
%   可选字段。如果某一条曲线想单独指定插值方式，可以加这个字段。
%   若不写，则 direct 使用 post.directThetaInterp，
%   WKB 使用 post.wkbThetaInterp。
%
% lineStyle, markerStyle:
%   可选字段。若不写，程序会自动分配不同线型和 marker。
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
        'label',      'WKB, N_\theta=17'), ...
    % struct( ...
    %     'dirName',    'direct_eps0p0001_dt0p0001_Nx40_Ntheta16_t2', ...
    %     'fileName',   'auto', ...
    %     'targetTime', 2, ...
    %     'solver',     'direct', ...
    %     'label',      'direct, N_\theta=16'), ...
    % struct( ...
    %     'dirName',    'direct_eps0p0001_dt0p0001_Nx40_Ntheta32_t2', ...
    %     'fileName',   'auto', ...
    %     'targetTime', 2, ...
    %     'solver',     'direct', ...
    %     'label',      'direct, N_\theta=32'), ...
    % struct( ...
    %     'dirName',    'direct_eps0p0001_dt0p0001_Nx40_Ntheta64_t2', ...
    %     'fileName',   'auto', ...
    %     'targetTime', 2, ...
    %     'solver',     'direct', ...
    %     'label',      'direct, N_\theta=64'), ...
    % struct( ...
    %     'dirName',    'direct_eps0p0001_dt0p0001_Nx40_Ntheta128_t2', ...
    %     'fileName',   'auto', ...
    %     'targetTime', 2, ...
    %     'solver',     'direct', ...
    %     'label',      'direct, N_\theta=128') ...
    };

%% ========================= 计算并画图 =========================

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [4, 4, 18, 12]);

ax = axes(fig);
hold(ax, 'on');
box(ax, 'on');
grid(ax, 'on');

fprintf('\n================ P(theta) 后处理 ================\n');
fprintf('direct 插值方式：%s\n', post.directThetaInterp);
fprintf('WKB 插值方式：   %s\n', post.wkbThetaInterp);
fprintf('细网格点数：     %d\n', post.NthetaPlot);
fprintf('每条曲线 marker 最多数：%d\n', post.maxRawMarkers);

for q = 1:numel(files)

    dataDir = fullfile(caseDir, files{q}.dirName);

    if strcmpi(files{q}.fileName, 'auto')
        if ~isfield(files{q}, 'targetTime') || isempty(files{q}.targetTime)
            error('fileName = auto 时必须设置 targetTime。');
        end

        filePath = local_find_snapshot_by_time( ...
            dataDir, files{q}.solver, files{q}.targetTime);

        fprintf('\n[自动选择] %s\n', filePath);
    else
        filePath = fullfile(dataDir, files{q}.fileName);

        if exist(filePath, 'file') ~= 2
            error('找不到文件：\n%s\n请检查 dirName 和 fileName 是否完整一致。', filePath);
        end
    end

    Sraw = load(filePath);
    S = local_select_data_struct(Sraw);

    method = local_get_reconstruction_method(files{q}, post);

    [thetaPlot, PPlot, thetaRaw, PRaw, info] = local_compute_Ptheta( ...
        S, files{q}.solver, method, post);

    [lineStyle, markerStyle] = local_plot_style(q, files{q});

    % 插值后的光滑曲线
    hLine = plot(ax, thetaPlot, PPlot, ...
        'LineStyle', lineStyle, ...
        'Marker', 'none', ...
        'LineWidth', post.lineWidth, ...
        'DisplayName', files{q}.label);

    % marker 画在插值前的原始点 thetaRaw, PRaw 上
    if post.showRawPoints
        markerIndex = local_select_raw_markers( ...
            thetaRaw, PRaw, post.maxRawMarkers, post.rawMarkerMode);

        if ~isempty(markerIndex)
            plot(ax, thetaRaw(markerIndex), PRaw(markerIndex), ...
                'LineStyle', 'none', ...
                'Marker', markerStyle, ...
                'MarkerSize', post.rawMarkerSize, ...
                'LineWidth', post.rawMarkerLineWidth, ...
                'Color', hLine.Color, ...
                'MarkerFaceColor', 'none', ...
                'HandleVisibility', 'off');
        end
    end

    fprintf('[Ptheta] %-24s solver=%-7s method=%-7s rawPeak=%.8f plotPeak=%.8f file=%s\n', ...
        files{q}.label, files{q}.solver, method, ...
        info.thetaPeakRaw, info.thetaPeakPlot, filePath);
end

xlabel(ax, '\theta', 'Interpreter', 'tex');

if post.normalizeP
    ylabel(ax, 'normalized P(\theta)', 'Interpreter', 'tex');
else
    ylabel(ax, 'P(\theta)', 'Interpreter', 'tex');
end

title(ax, 'Trait marginal comparison', 'Interpreter', 'tex');
legend(ax, 'Location', 'best', 'Box', 'off');
xlim(ax, [0, 1]);

set(ax, 'FontSize', 12);

pngFile = fullfile(post.outputDir, [post.figureName, '.png']);
figFile = fullfile(post.outputDir, [post.figureName, '.fig']);

try
    exportgraphics(fig, pngFile, 'Resolution', 300);
catch
    saveas(fig, pngFile);
end

savefig(fig, figFile);

fprintf('\n图片已保存：\n%s\n%s\n', pngFile, figFile);

%% ========================================================================
% 局部函数
% ========================================================================

function method = local_get_reconstruction_method(oneFile, post)
% 根据 solver 类型选择插值方式。
% 若 files 中单独给了 reconMethod，则优先使用该方法。

if isfield(oneFile, 'reconMethod') && ~isempty(oneFile.reconMethod)
    method = oneFile.reconMethod;
    return;
end

switch lower(oneFile.solver)
    case 'direct'
        method = post.directThetaInterp;
    case 'wkb'
        method = post.wkbThetaInterp;
    otherwise
        error('未知 solver 类型：%s。', oneFile.solver);
end
end

function filePath = local_find_snapshot_by_time(dataDir, solver, targetTime)
% 在 dataDir 中查找最接近 targetTime 的 snapshot 文件。
% 兼容：
%   snapshot_direct_t2.mat
%   snapshot_direct_t0p5.mat
%   snapshot_direct_0002000_t2.mat
%   snapshot_direct_step0002000_t2.mat

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
%   snapshot_direct_t2.mat
%   snapshot_direct_t0p5.mat
%   snapshot_direct_t2p5em4.mat
%   snapshot_direct_0002000_t2.mat
%   snapshot_direct_step0002000_t2.mat

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
% 兼容不同保存格式：
%   snapshot
%   direct
%   wkb
%   one
%   或者文件中只有一个结构体变量

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

function [thetaPlot, PPlot, thetaRaw, PRaw, info] = local_compute_Ptheta(S, solver, method, post)
% 根据 solver 类型计算 P(theta)。
%
% direct:
%   nRaw -> 插值 nFine -> PPlot
%
% WKB:
%   WRaw, uRaw -> 分别插值 WFine, UFine
%   nFine = WFine exp(UFine/eps) -> PPlot

[x, thetaRaw] = local_get_grid(S);

Nx = numel(x);
Ntheta = numel(thetaRaw);

thetaPlot = linspace(0, 1, post.NthetaPlot + 1).';
thetaPlot(end) = [];

switch lower(solver)

    case 'direct'
        if ~isfield(S, 'n') || isempty(S.n)
            error('direct 文件中没有变量 n。');
        end

        nRaw = local_to_xt_matrix(S.n, Nx, Ntheta, 'n');

        PRaw = trapz(x(:), nRaw, 1).';
        PRaw = real(PRaw(:));

        if strcmpi(method, 'none')
            thetaPlot = thetaRaw;
            PPlot = PRaw;
        else
            nFine = local_reconstruct_theta_matrix( ...
                thetaRaw, nRaw, thetaPlot, method);

            nFine = real(nFine);

            if post.clipNegative
                nFine(nFine < 0) = 0;
            end

            PPlot = trapz(x(:), nFine, 1).';
            PPlot = real(PPlot(:));
        end

    case 'wkb'
        if ~isfield(S, 'W') || isempty(S.W) || ~isfield(S, 'u') || isempty(S.u)
            error('WKB 文件中需要同时包含 W 和 u。');
        end

        WRaw = local_to_xt_matrix(S.W, Nx, Ntheta, 'W');
        uRaw = S.u;

        if isvector(uRaw)
            URaw = repmat(real(uRaw(:).'), Nx, 1);
        else
            URaw = local_to_xt_matrix(uRaw, Nx, Ntheta, 'u');
        end

        epsVal = local_get_eps(S);

        nRaw = max(real(WRaw), 0) .* exp(min(max(real(URaw) / epsVal, -700), 700));
        PRaw = trapz(x(:), nRaw, 1).';
        PRaw = real(PRaw(:));

        if strcmpi(method, 'none')
            thetaPlot = thetaRaw;
            PPlot = PRaw;
        else
            WFine = local_reconstruct_theta_matrix( ...
                thetaRaw, WRaw, thetaPlot, method);

            UFine = local_reconstruct_theta_matrix( ...
                thetaRaw, URaw, thetaPlot, method);

            WFine = real(WFine);
            UFine = real(UFine);

            if post.clipNegative
                WFine(WFine < 0) = 0;
            end

            expo = UFine / epsVal;
            expo = min(max(expo, -700), 700);

            nFine = WFine .* exp(expo);

            PPlot = trapz(x(:), nFine, 1).';
            PPlot = real(PPlot(:));
        end

    otherwise
        error('未知 solver 类型：%s。只能是 direct 或 wkb。', solver);
end

if post.normalizeP
    PRaw = local_normalize_periodic(thetaRaw, PRaw);
    PPlot = local_normalize_periodic(thetaPlot, PPlot);
end

[~, idRaw] = max(PRaw);
[~, idPlot] = max(PPlot);

info.thetaPeakRaw = thetaRaw(idRaw);
info.thetaPeakPlot = thetaPlot(idPlot);
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

function AFine = local_reconstruct_theta_matrix(theta, A, thetaFine, method)
% 对矩阵 A(x,theta) 的每一行在 theta 方向做周期重构。

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
% 周期重构接口。

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
% interp1 周期插值。
% 使用左右各补一个周期点，避免 thetaRaw 不从 0 开始时出现 NaN。

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
% 五阶 WENO 型周期插值。
%
% 说明：
%   这是后处理画图用的 WENO 型重构，不改变核心求解器。
%   若网格太少，则自动退化为 pchip。

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
% 周期下标，返回 1,...,K。
id = mod(i - 1, K) + 1;
end

function P = local_normalize_periodic(theta, P)
% 周期区间上的归一化，使 int_0^1 P(theta) dtheta = 1。

theta = theta(:);
P = real(P(:));

if numel(theta) > 1
    mass = trapz([theta; theta(1) + 1], [P; P(1)]);
else
    mass = sum(P);
end

if isfinite(mass) && abs(mass) > realmin
    P = P / mass;
end
end

function [lineStyle, markerStyle] = local_plot_style(q, oneFile)
% 为第 q 条曲线选择线型和 marker。
% 若 files 中单独设置 lineStyle 或 markerStyle，则优先使用。

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

function markerIndex = local_select_raw_markers(thetaRaw, PRaw, maxRawMarkers, mode)
% 从插值前的原始 theta 网格中选择 marker 位置。
%
% thetaRaw, PRaw:
%   插值前的原始数据。
%
% maxRawMarkers:
%   每条曲线最多显示几个 marker。
%
% mode:
%   'uniform'
%       均匀抽样原始点。
%
%   'uniform_peak'
%       均匀抽样，同时尽量包含峰值点。
%
%   'all'
%       显示全部原始点，不受 maxRawMarkers 限制。

thetaRaw = thetaRaw(:);
PRaw = PRaw(:);

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

        [~, peakIndex] = max(PRaw);

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