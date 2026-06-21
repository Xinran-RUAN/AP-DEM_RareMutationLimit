function H = solve_H_fd_1d(D, Kx, rho, op)
%SOLVE_H_FD_1D 计算每个 theta_k 对应的离散主特征值 H_k。
%
% 离散特征值问题为
%   -D_k L_x N - diag(K-rho) N = H_k N.
%
% 速度说明：
%   H 的计算通常是 WKB 主程序里最耗时的部分。每一步都要对所有 theta_k
%   求一次空间主特征值。因此，Nx 和 Ntheta 同时加大时，运行时间会明显增加。
%   这里做了两个优化：
%     1. 预先构造梯形权重下的对称相似变换后的 Lx，避免在 theta 循环中重复做
%        S*A*S^{-1}；
%     2. 对 eigs 使用上一个 theta 的特征向量作为 v0 初值，因为相邻 theta 的
%        D(theta) 接近，主特征向量通常也接近；小矩阵时直接用 dense eig，避免
%        eigs 的迭代启动开销。
%
% 注意：这里没有改变离散模型，只是加速相同离散特征值问题的求解。

Ktheta = numel(D);
nx = op.nx;
H = zeros(1, Ktheta);

% Neumann reflected Laplacian 在梯形权重内积下自伴。用相似变换转成欧氏
% 内积下的对称矩阵，随后调用 eig/eigs 的最小代数特征值。
V = Kx(:) - rho(:);
wx = utils.trapz_weights(nx, op.dx);
S = spdiags(sqrt(wx), 0, nx, nx);
Sinv = spdiags(1 ./ sqrt(wx), 0, nx, nx);
Lsym = S * op.Lx * Sinv;
Lsym = (Lsym + Lsym')/2;
Vdiag = spdiags(V, 0, nx, nx);

% 小矩阵用 dense eig 反而更快；较大矩阵用 eigs。
denseThreshold = 150;
useDense = (nx <= denseThreshold);

opts = struct();
opts.tol = 1e-10;
opts.maxit = 300;
v0 = [];

for k = 1:Ktheta
    As = -D(k) * Lsym - Vdiag;
    As = (As + As')/2;  % 消除舍入误差导致的非对称小扰动

    if useDense
        ev = eig(full(As));
        H(k) = min(real(ev));
        continue;
    end

    try
        if ~isempty(v0)
            opts.v0 = v0;
        elseif isfield(opts, 'v0')
            opts = rmfield(opts, 'v0');
        end
        [v, lam] = eigs(As, 1, 'sa', opts);
        H(k) = real(lam(1,1));
        v0 = v;
    catch
        try
            if ~isempty(v0)
                opts.v0 = v0;
            elseif isfield(opts, 'v0')
                opts = rmfield(opts, 'v0');
            end
            [v, lam] = eigs(As, 1, 'smallestreal', opts);
            H(k) = real(lam(1,1));
            v0 = v;
        catch
            ev = eig(full(As));
            H(k) = min(real(ev));
            v0 = [];
        end
    end
end
end
