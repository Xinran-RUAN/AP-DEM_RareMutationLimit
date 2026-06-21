function [Wnew, info] = update_w_full_eigen_relax_fd_1d(W, u, D, Kx, rho, Hupdate, Hraw, eps, dt, op, par)
%UPDATE_W_FULL_EIGEN_RELAX_FD_1D 双网格 full-eigen-relax 的 W 更新。
%
%  重要：本函数没有把 W 方程近似成逐 theta 的 x-特征向量投影。
%  实际求解的仍然是完整的 WKB 振幅一步格式：
%
%    eps (W^{n+1}-W^n)/dt
%      - D(theta) delta_x^2 W^{n+1}
%      - eps^2 delta_theta^2 W^{n+1}
%      + 2 eps D_theta F(W^{n+1},u^{n+1})
%    = W^{n+1} (K-rho^n+H^n-2 eps delta_theta^2 u^{n+1}).
%
%  因此 x 扩散、theta 扩散和 theta 输运项是在同一个线性系统中同时处理的。
%  与旧 split 更新的区别只在于：full-eigen-relax 分支要求 Hraw 最好由 W
%  粗 theta 网格上的同一个 x-算子直接计算得到。这样当 eps 很小时，完整格式
%  的主导 x-方向快松弛自然具有
%
%       (-D_k Lx - diag(K-rho) - Hraw_k I) N_k = 0
%
%  的特征向量结构。
%
%  输入：
%    W       : nx x Ktheta_W 的粗网格振幅旧值
%    u       : 1  x Ktheta_W 的粗网格相位新值，通常由细 u 插值得到
%    D       : 1  x Ktheta_W 的扩散率
%    Kx,rho  : nx x 1
%    Hupdate : W 方程中实际使用的 H。通常取 Hraw；若需要和 H gauge 完全一致，
%              可在外部传入 shifted H。
%    Hraw    : 未做 Hamiltonian gauge shift 的主特征值，只用于诊断特征向量残差。
%    eps,dt,op,par : 参数和算子
%
%  输出：
%    Wnew : nx x Ktheta_W
%    info : 诊断信息，包括 xEigenResidualInfMax。该残差不是格式残差，而是
%           检查 Wnew 每个 theta 截面是否接近对应 x 主特征向量方向。

if nargin < 11 || isempty(par)
    par = struct();
end
if nargin < 7 || isempty(Hraw)
    Hraw = Hupdate;
end

info = struct();
info.amplitudeVariantUsed = 'full-eigen-relax';
info.fullWUpdateKeepsThetaTerms = true;
info.comment = 'full W update: x diffusion + theta diffusion + theta transport in one linear system';

Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;

D = real(D(:).');
Kx = real(Kx(:));
rho = real(rho(:));
Hupdate = real(Hupdate(:).');
Hraw = real(Hraw(:).');
u = real(u(:).');

if size(W,1) ~= nx || size(W,2) ~= Ktheta
    error('update_w_full_eigen_relax_fd_1d:BadWSize', ...
        'W 的尺寸应为 op.nx x op.Ntheta = %d x %d。', nx, Ktheta);
end
if numel(D) ~= Ktheta || numel(Hupdate) ~= Ktheta || numel(Hraw) ~= Ktheta || numel(u) ~= Ktheta
    error('update_w_full_eigen_relax_fd_1d:BadThetaSize', ...
        'D/H/u 的长度必须等于 op.Ntheta。');
end
if numel(Kx) ~= nx || numel(rho) ~= nx
    error('update_w_full_eigen_relax_fd_1d:BadXSize', ...
        'Kx/rho 的长度必须等于 op.nx。');
end

%-------------------------
% 完整 W 方程矩阵。
% 与 update_w_split_fd_1d 的离散形式保持一致，只是这里额外输出诊断，
% 并明确允许 Hupdate 来自 W 粗网格的 raw principal eigenvalue。
%-------------------------
[Ifix, DblockLx, Lth] = local_fixed_matrices(D, op);

ddu = (op.Ltheta * u(:)).';
extra = -2*eps*ddu;
r = src.reaction_vector_1d(Kx, rho, Hupdate, extra, op);

A = eps*Ifix - dt*DblockLx - dt*eps^2*Lth - dt*spdiags(r,0,N,N);

T = src.upwind_transport_matrix(u, op);
if ~isfield(par,'transportImplicit') || par.transportImplicit
    A = A + dt*2*eps*T;
    info.transportImplicitUsed = true;
    b = eps*W(:);
else
    b = eps*W(:) - dt*2*eps*(T*W(:));
    info.transportImplicitUsed = false;
end

Wvec = A \ b;
Wnew = reshape(real(Wvec), nx, Ktheta);

info.numNonfiniteW = nnz(~isfinite(Wnew(:)));
if info.numNonfiniteW > 0
    Wnew(~isfinite(Wnew)) = 0;
end

info.minWBeforeClip = min(Wnew(:));
info.maxWBeforeClip = max(Wnew(:));
if local_get_bool(par, 'clipNegativeW', true)
    floorW = local_get_num(par, 'WFloor', 0);
    info.numNegativeW = nnz(Wnew < floorW);
    Wnew = max(Wnew, floorW);
else
    info.numNegativeW = nnz(Wnew < 0);
end

info.WmaxAfter = max(Wnew(:));
info.WminAfter = min(Wnew(:));
info.matrixDiagMin = min(abs(diag(A)));
info.matrixDiagMax = max(abs(diag(A)));
info.matrixCondEst = NaN;
info.matrixRcondEst = NaN;
if local_get_bool(par, 'monitorMatrixCondition', false)
    try
        info.matrixCondEst = condest(A);
        info.matrixRcondEst = 1 / info.matrixCondEst;
    catch
        info.matrixCondEst = NaN;
        info.matrixRcondEst = NaN;
    end
end

%-------------------------
% 线性系统残差。
%-------------------------
rhsNorm = max(norm(b, inf), realmin);
info.solveResidualInf = norm(A*Wnew(:) - b, inf) / rhsNorm;

%-------------------------
% 小 eps 特征向量流形诊断。
% 对每个 theta_k，检查 x 方向主导算子残差：
%   B_k W_k = (-D_k Lx - diag(K-rho) - Hraw_k I) W_k.
% 该诊断不参与格式，只用于判断 W 是否已经接近对应主特征向量方向。
%-------------------------
Ix = op.Ix;
Vdiag = spdiags(Kx - rho, 0, nx, nx);
xEigRes = zeros(1, Ktheta);
xEigResMass = zeros(1, Ktheta);
for k = 1:Ktheta
    Bk = -D(k)*op.Lx - Vdiag - Hraw(k)*Ix;
    wk = Wnew(:,k);
    xEigRes(k) = norm(Bk*wk, inf) / max(norm(wk, inf), realmin);
    xEigResMass(k) = norm(Bk*wk, 2) / max(norm(wk, 2), realmin);
end
info.xEigenResidualInf = xEigRes;
info.xEigenResidualInfMax = max(xEigRes);
info.xEigenResidualL2Max = max(xEigResMass);
info.HupdateMinusHrawMaxAbs = max(abs(Hupdate(:)-Hraw(:)));
info.failed = false;
info.failureReason = '';
if info.numNonfiniteW > 0
    info.failed = true;
    info.failureReason = 'W 线性求解产生 NaN 或 Inf，已置零后继续。';
end
end

function [I, DblockLx, Lth] = local_fixed_matrices(D, op)
% Persistent cache keyed by grid size and D vector.
persistent cache

D = real(D(:));
key = sprintf('nx%d_K%d_dx%.17g_dt%.17g_D%.17g_%.17g_%.17g', ...
    op.nx, op.Ntheta, op.dx, op.dtheta, sum(D), sum(D.^2), sum(abs(D)));

if isempty(cache) || ~isfield(cache, 'key') || ~strcmp(cache.key, key)
    Ktheta = op.Ntheta;
    nx = op.nx;
    cache = struct();
    cache.key = key;
    cache.I = speye(nx*Ktheta);
    cache.DblockLx = kron(spdiags(D,0,Ktheta,Ktheta), op.Lx);
    cache.Lth = op.Ltheta_big;
end

I = cache.I;
DblockLx = cache.DblockLx;
Lth = cache.Lth;
end

function tf = local_get_bool(s, name, defaultValue)
tf = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end

function val = local_get_num(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end
