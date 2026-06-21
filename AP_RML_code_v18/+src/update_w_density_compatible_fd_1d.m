function [Wnew, info] = update_w_density_compatible_fd_1d(W, u, D, Kx, rho, H, eps, dt, op, par)
%UPDATE_W_DENSITY_COMPATIBLE_FD_1D 旧的半离散 density-compatible 全离散振幅更新。
%   注意：小 epsilon 正式计算默认不使用该矩阵版本；它仅保留作对比。该版本直接离散 (E^{-1} delta_theta^2(W E))，其中 E=exp(u/eps)。
%   因而重构密度 n=W*E 在离散层面严格满足 conservative density balance，
%   更适合验证文章中 weighted remainder 的理论条件。
info = struct();
Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;

DblockLx = kron(spdiags(D(:),0,Ktheta,Ktheta), op.Lx);
Btheta = src.density_compatible_theta_matrix(u, eps, op);

ddu = (op.Ltheta * u(:)).';
Hhat = src.numerical_hamiltonian(u, op.dtheta, par.phaseHamiltonian, par.lfAlpha);
extra = Hhat - eps*ddu;
r = src.reaction_vector_1d(Kx, rho, H, extra, op);

A = eps*speye(N) - dt*DblockLx - dt*eps^2*Btheta - dt*spdiags(r,0,N,N);
b = eps*W(:);
Wnew = reshape(A \ b, nx, Ktheta);
Wnew = real(Wnew);
info.WmaxAfter = max(Wnew(:));
info.WminAfter = min(Wnew(:));
end
