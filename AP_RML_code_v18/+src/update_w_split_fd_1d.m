function [Wnew, info] = update_w_split_fd_1d(W, u, D, Kx, rho, H, eps, dt, op, par)
%UPDATE_W_SPLIT_FD_1D split-WKB 振幅方程的一步更新。
%   默认格式：反应项、x/theta 扩散项隐式；theta 输运项可隐式或显式。
%   这里组装的是关于 W^{n+1} 的线性系统。
%
% 速度优化：固定的 D-block Lx、theta Laplacian、单位矩阵用 persistent cache，
% 每步只更新 reaction diagonal 和 u 依赖的 upwind transport 矩阵。

info = struct();
Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;

[Ifix, DblockLx, Lth] = local_fixed_matrices(D, op);

ddu = (op.Ltheta * u(:)).';
extra = -2*eps*ddu;
r = src.reaction_vector_1d(Kx, rho, H, extra, op);

A = eps*Ifix - dt*DblockLx - dt*eps^2*Lth - dt*spdiags(r,0,N,N);

T = src.upwind_transport_matrix(u, op);
if ~isfield(par,'transportImplicit') || par.transportImplicit
    A = A + dt*2*eps*T;
    b = eps*W(:);
else
    b = eps*W(:) - dt*2*eps*(T*W(:));
end

Wnew = reshape(A \ b, nx, Ktheta);
Wnew = real(Wnew);
info.WmaxAfter = max(Wnew(:));
info.WminAfter = min(Wnew(:));
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
