function some_try
% --------- 参数设置 ----------
tic 
N = 10000;                   % 网格数（可改为更大）
h = 1 / (N + 1);            % 网格步长
x = linspace(0, 1, N+2);    % 包含边界点
e = ones(N,1);

% --------- 构造有限差分矩阵（Dirichlet 边界） ----------
A = spdiags([-e 2*e -e], -1:1, N, N) / h^2;

% --------- Shift-and-Invert 反幂法求最小特征值 ----------
sigma = 0;                          % shift（目标最小值附近）
S = A - sigma * speye(N);           % shifted 矩阵
[L,U] = lu(S);                      % LU 分解（一次性）

w = rand(N,1);                      % 初始向量
w = w / norm(w);
max_iter = 100;
tol = 1e-10;

for k = 1:max_iter
    w_old = w;
    w = U \ (L \ w);               % 解 (A - sigma*I)x = b
    w = w / norm(w);               % 单位化
    if norm(w - w_old) < tol
        break;
    end
end

lambda = (w' * A * w) / (w' * w);  % Rayleigh 商估计特征值

% --------- 输出结果 ----------
fprintf("最小特征值估计: lambda = %.12f\n", lambda);
fprintf("迭代次数: %d\n", k);
toc
end

