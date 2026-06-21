function op = build_operators_1d(grid)
%BUILD_OPERATORS_1D Sparse finite-difference operators.
% x-grid: endpoints x_0,...,x_N with reflected Neumann ghost values.
% theta-grid: periodic nodes theta_k = k/K.

nx = grid.Nx + 1;
K = grid.Ntheta;
dx = grid.dx;
dth = grid.dtheta;

ex = ones(nx,1);
Lx = spdiags([ex -2*ex ex], [-1 0 1], nx, nx) / dx^2;
% Reflected endpoint Neumann: u_{-1}=u_1 and u_{N+1}=u_{N-1}.
Lx(1,:) = 0;     Lx(1,1) = -2/dx^2; Lx(1,2) =  2/dx^2;
Lx(end,:) = 0;   Lx(end,end) = -2/dx^2; Lx(end,end-1) = 2/dx^2;

et = ones(K,1);
Ltheta = spdiags([et -2*et et], [-1 0 1], K, K) / dth^2;
Ltheta(1,K) = 1/dth^2;
Ltheta(K,1) = 1/dth^2;

op.nx = nx;
op.Ntheta = K;
op.x = grid.x;
op.theta = grid.theta;
op.dx = dx;
op.dtheta = dth;
op.Lx = sparse(Lx);
op.Ltheta = sparse(Ltheta);
op.Ix = speye(nx);
op.Itheta = speye(K);
op.Ltheta_big = kron(op.Ltheta, op.Ix);
end
