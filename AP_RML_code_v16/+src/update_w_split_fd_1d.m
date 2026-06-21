function [Wnew, info] = update_w_split_fd_1d(W, u, D, Kx, rho, H, eps, dt, op, par)
%UPDATE_W_SPLIT_FD_1D split-WKB 振幅方程的一步更新。
%   默认格式：反应项、x/theta 扩散项隐式；theta 输运项可隐式或显式。
%   这里组装的是关于 W^{n+1} 的线性系统。
info = struct();
Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;

% x 方向扩散：每个 theta_k 上的系数为 D(theta_k)。
DblockLx = kron(spdiags(D(:),0,Ktheta,Ktheta), op.Lx);

% theta 方向扩散和相位曲率修正。
Lth = op.Ltheta_big;
ddu = (op.Ltheta * u(:)).';
extra = -2*eps*ddu;
r = src.reaction_vector_1d(Kx, rho, H, extra, op);

% 线性系统 A W^{n+1} = b。
A = eps*speye(N) - dt*DblockLx - dt*eps^2*Lth - dt*spdiags(r,0,N,N);

% 保守 upwind 离散 D_theta F, F=-W u_theta。
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
