function [H, info] = solve_H_fd_1d(D, Kx, rho, op, par)
%SOLVE_H_FD_1D 计算每个 theta_k 对应的离散主特征值 H_k。
%
% 离散特征值问题为
%   -D_k L_x N - diag(K-rho) N = H_k N.
%
% 新增高效模式（用户接口字段 par.HSolverMode）：
%   'direct-all-theta'          每个 theta 都直接求特征值；
%   'direct-unique-D'           对 D(theta) 去重后求特征值，默认推荐；
%   'D-pchip'                   固定 rho 后，把 H 看成 H(D;rho)，在少量 D 点求解后 pchip 插值；
%   'D-pchip-local-correct'     D-pchip 后，对 par.HLocalCorrectionIndices 指定的峰附近点直接校正。
%
% 其它字段：
%   par.HDtol       D 去重容差，默认 1e-12；
%   par.HInterpND   D-pchip 采样点数，默认 min(24,numel(D))；
%   par.eigenSolver {'auto','eig-full','eigs'}。
%
% 保持旧调用兼容：H = src.solve_H_fd_1d(D,Kx,rho,op)。

if nargin < 5 || isempty(par)
    par = struct();
end

D = real(D(:));
Ktheta = numel(D);
nx = op.nx;
H = zeros(1, Ktheta);

mode = local_get_string(par, 'HSolverMode', 'direct-unique-D');
mode = lower(strrep(mode, '_', '-'));

info = struct();
info.mode = mode;
info.numTheta = Ktheta;
info.numEigenSolves = 0;
info.numDirectCorrections = 0;
info.usedInterpolation = false;

ctx = local_make_eigen_context(Kx, rho, op, par);

switch mode
    case {'direct-all-theta','all-theta','direct'}
        for k = 1:Ktheta
            H(k) = local_min_eigen_for_D(D(k), ctx);
        end
        info.numEigenSolves = Ktheta;

    case {'direct-unique-d','unique-d','unique'}
        tolD = local_get_field(par, 'HDtol', 1e-12);
        if tolD <= 0 || ~isfinite(tolD)
            tolD = 1e-12;
        end
        Dkey = round(D / tolD) * tolD;
        [DuniqKey, ~, ic] = unique(Dkey, 'stable'); %#ok<ASGLU>
        % 用每个 key 第一次出现的原始 D 值，避免 round 后轻微改动数值。
        Duniq = zeros(numel(DuniqKey), 1);
        for r = 1:numel(DuniqKey)
            first = find(ic == r, 1, 'first');
            Duniq(r) = D(first);
        end
        Huniq = zeros(numel(Duniq), 1);
        for r = 1:numel(Duniq)
            Huniq(r) = local_min_eigen_for_D(Duniq(r), ctx);
        end
        H = Huniq(ic).';
        info.numEigenSolves = numel(Duniq);
        info.numUniqueD = numel(Duniq);

    case {'d-pchip','dpchip','dinterp','d-interp','d-pchip-local-correct','dinterp-local-correct'}
        ND = round(local_get_field(par, 'HInterpND', min(24, Ktheta)));
        ND = max(2, min(ND, max(2, Ktheta)));
        Dmin = min(D);
        Dmax = max(D);

        if abs(Dmax - Dmin) <= eps(max(1, abs(Dmin)))
            H(:) = local_min_eigen_for_D(Dmin, ctx);
            info.numEigenSolves = 1;
            info.usedInterpolation = false;
        else
            Dsample = linspace(Dmin, Dmax, ND).';
            Hsample = zeros(ND, 1);
            for r = 1:ND
                Hsample(r) = local_min_eigen_for_D(Dsample(r), ctx);
            end
            H = interp1(Dsample, Hsample, D, 'pchip', 'extrap').';
            info.numEigenSolves = ND;
            info.usedInterpolation = true;
            info.HInterpND = ND;

            if any(strcmp(mode, {'d-pchip-local-correct','dinterp-local-correct'}))
                idx = [];
                if isfield(par, 'HLocalCorrectionIndices') && ~isempty(par.HLocalCorrectionIndices)
                    idx = unique(round(par.HLocalCorrectionIndices(:).'));
                    idx = idx(idx >= 1 & idx <= Ktheta);
                end
                if ~isempty(idx)
                    for kk = idx
                        H(kk) = local_min_eigen_for_D(D(kk), ctx);
                    end
                    info.numEigenSolves = info.numEigenSolves + numel(idx);
                    info.numDirectCorrections = numel(idx);
                    info.localCorrectionIndices = idx;
                end
            end
        end

    otherwise
        error('solve_H_fd_1d:UnknownMode', '未知 HSolverMode: %s', mode);
end

H = real(H(:).');
end

% =====================================================================
% 局部函数
% =====================================================================

function ctx = local_make_eigen_context(Kx, rho, op, par)
nx = op.nx;
V = Kx(:) - rho(:);
wx = utils.trapz_weights(nx, op.dx);
S = spdiags(sqrt(wx), 0, nx, nx);
Sinv = spdiags(1 ./ sqrt(wx), 0, nx, nx);
Lsym = S * op.Lx * Sinv;
Lsym = (Lsym + Lsym')/2;
Vdiag = spdiags(V, 0, nx, nx);

ctx = struct();
ctx.nx = nx;
ctx.Lsym = Lsym;
ctx.Vdiag = Vdiag;
ctx.eigenSolver = lower(local_get_string(par, 'eigenSolver', 'auto'));
ctx.denseThreshold = local_get_field(par, 'eigenDenseThreshold', 180);
ctx.opts = struct('tol', local_get_field(par, 'eigenTol', 1e-10), ...
                  'maxit', local_get_field(par, 'eigenMaxit', 300), ...
                  'isreal', true, 'issym', true);
end

function lam = local_min_eigen_for_D(Dval, ctx)
As = -Dval * ctx.Lsym - ctx.Vdiag;
As = (As + As')/2;

solver = ctx.eigenSolver;
if strcmp(solver, 'auto')
    if ctx.nx <= ctx.denseThreshold
        solver = 'eig-full';
    else
        solver = 'eigs';
    end
end

switch solver
    case {'eig-full','eig','dense'}
        ev = eig(full(As));
        lam = min(real(ev));

    case {'eigs','sparse'}
        try
            lam = eigs(As, 1, 'sa', ctx.opts);
            lam = real(lam(1,1));
        catch
            try
                lam = eigs(As, 1, 'smallestreal', ctx.opts);
                lam = real(lam(1,1));
            catch
                ev = eig(full(As));
                lam = min(real(ev));
            end
        end

    otherwise
        error('solve_H_fd_1d:UnknownEigenSolver', '未知 eigenSolver: %s', ctx.eigenSolver);
end
end

function v = local_get_field(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    v = double(s.(name));
end
end

function v = local_get_string(s, name, defaultValue)
v = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = char(s.(name));
end
end
