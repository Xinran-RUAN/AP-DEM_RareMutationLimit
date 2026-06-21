function [W, u, D, K, grid] = initial_wkb_1d(par)
%INITIAL_WKB_1D WKB 求解器的初值和系数。
%
% 初值直接由主程序中的函数句柄给出：
%
%     par.u0_fun(theta)
%     par.W0_fun(x,theta)
%
% 其中 theta 是 1 x Ntheta 行向量，x 是 (Nx+1) x 1 列向量。
% 推荐在主程序中直接写表达式，例如
%
%     par.u0_fun = @(theta) -0.05*(1-cos(2*pi*(theta-0.70)));
%     par.W0_fun = @(x,theta) 1 + 0.*x + 0.*theta;

grid = model.make_grid_1d(par);
x = grid.x;
theta = grid.theta;

D = par.D_fun(theta);
K = par.K_fun(x);

%==========================================================================
% 1. 初值
%==========================================================================

if ~isfield(par, 'u0_fun') || isempty(par.u0_fun)
    u = zeros(1, par.Ntheta);
else
    u = par.u0_fun(theta);
end

if ~isfield(par, 'W0_fun') || isempty(par.W0_fun)
    W = ones(numel(x), par.Ntheta);
else
    W = par.W0_fun(x, theta);
end

%==========================================================================
% 2. 整理尺寸
%==========================================================================

u = real(u);

if isscalar(u)
    u = u * ones(1, par.Ntheta);
end

u = u(:).';

if numel(u) ~= par.Ntheta
    error('u0_fun 返回的长度为 %d，但 par.Ntheta=%d。', numel(u), par.Ntheta);
end

W = real(W);

if isscalar(W)
    W = W * ones(numel(x), par.Ntheta);
elseif isvector(W)
    Wv = W(:);

    if numel(Wv) == numel(x)
        W = repmat(Wv, 1, par.Ntheta);
    elseif numel(Wv) == par.Ntheta
        W = repmat(Wv(:).', numel(x), 1);
    else
        error('W0_fun 返回了长度为 %d 的向量，无法匹配 x 或 theta 网格。', numel(Wv));
    end
elseif isequal(size(W), [par.Ntheta, numel(x)])
    W = W.';
end

if ~isequal(size(W), [numel(x), par.Ntheta])
    error('W0_fun 返回的尺寸为 %s，但应为 [%d,%d]。', ...
        mat2str(size(W)), numel(x), par.Ntheta);
end

if any(~isfinite(u(:)))
    error('初始相位 u 含有 NaN 或 Inf。');
end

if any(~isfinite(W(:)))
    error('初始振幅 W 含有 NaN 或 Inf。');
end

if any(W(:) < 0)
    error('初始振幅 W 必须非负。');
end

%==========================================================================
% 3. 初始 gauge
%==========================================================================

gaugeMode = 'max';
if isfield(par, 'initialGaugeMode') && ~isempty(par.initialGaugeMode)
    gaugeMode = char(par.initialGaugeMode);
end

switch lower(gaugeMode)

    case {'max','legacy-max','max-zero'}
        [W, u] = src.normalize_gauge(W, u, par.eps, 'max');

    case {'none','off'}
        % 不做 gauge。

    otherwise
        error('未知 initialGaugeMode：%s。', gaugeMode);
end

%==========================================================================
% 4. 可选初始质量缩放
%==========================================================================

if isfield(par, 'prepareInitialMass') && logical(par.prepareInitialMass)

    op = src.build_operators_1d(grid);

    rhoMode = 'direct-log';
    if isfield(par, 'rhoReconstruction') && ~isempty(par.rhoReconstruction)
        rhoMode = char(par.rhoReconstruction);
    end

    [rhoCurrent, ~] = src.reconstruct_rho_1d(W, u, par.eps, op, rhoMode);

    rhoTarget = max(0.5 * K(:), 1e-12);
    fac = rhoTarget(:) ./ max(rhoCurrent(:), realmin);

    W = bsxfun(@times, W, fac);
end

end