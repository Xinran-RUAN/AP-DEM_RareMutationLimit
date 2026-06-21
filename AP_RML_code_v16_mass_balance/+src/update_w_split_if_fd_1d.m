function [Wnew, info] = update_w_split_if_fd_1d(W, uOld, uNew, D, Kx, rho, H, eps, dt, op, par)
%UPDATE_W_SPLIT_IF_FD_1D 纯 WKB 框架下的时间平衡 split 振幅更新。
%
%  这个函数仍然求解振幅变量 W，而不是把问题改写成原始密度 n 的求解。
%  它的目的，是去掉旧 split 全离散中最刚性的时间链式法则误差：
%
%       n = W exp(u/eps).
%
%  旧格式使用 eps*(W^{n+1}-W^n)/dt。这里改为
%
%       eps*(W^{n+1}-q_k W^n)/dt,
%       q_k = exp((u^n_k-u^{n+1}_k)/eps),
%  这相当于在 WKB 变量内部加入时间方向 integrating factor。注意这里
%  只使用同一个 theta 点上的 q_k，不出现
%       exp((u_{k+1}-u_k)/eps)
%  这种相邻点指数比值，因此不会产生 density-compatible-if 的奇异矩阵。
%
%  经过相位方程消去 H 项后，振幅方程中的源项写成
%       r_if = K-rho - Hhat(u^n) - eps*delta_{theta theta} u^{n+1},
%  其中 Hhat 近似 -|u_theta|^2。反应项默认采用 Patankar 型处理，
%  以避免正增长项在小 eps 时破坏非负性。

if nargin < 11 || isempty(par)
    par = struct();
end

Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;
info = struct();

reactionMode = local_get_string(par, 'reactionDiscretization', 'patankar');
expClip = local_get_num(par, 'ifExpClip', local_get_num(par, 'expClip', 80));

%--------------------------------------------------------------
% 1. q_k = exp((uOld-uNew)/eps)。这是对时间方向 WKB 链式法则的补偿。
%    为了防止浮点溢出，做指数裁剪；裁剪次数会写入 info 供诊断。
%--------------------------------------------------------------
qLog = (uOld(:).' - uNew(:).') / eps;
qLogClipped = min(max(qLog, -expClip), expClip);
q = exp(qLogClipped);
qW = bsxfun(@times, W, q);

info.qLogMin = min(qLog);
info.qLogMax = max(qLog);
info.qLogAbsMax = max(abs(qLog));
info.qClipCount = sum(qLog ~= qLogClipped);
info.qClipFraction = info.qClipCount / max(numel(qLog),1);

%--------------------------------------------------------------
% 2. 组装 W 变量上的 split-WKB 线性系统。
%    注意这里没有 E_{k+1}/E_k，因此不会在 eps 很小时因相邻相位差产生奇异矩阵。
%--------------------------------------------------------------
DblockLx = kron(spdiags(D(:),0,Ktheta,Ktheta), op.Lx);
Lth = op.Ltheta_big;
T = src.upwind_transport_matrix(uNew, op);

% Hhat 使用相位方程中的显式 Hamiltonian，即旧相位 uOld。
% 常数 gauge 不影响差分，因此这里用 uOld 或同 gauge 的 uOld 均可。
HhatOld = src.numerical_hamiltonian(uOld, op.dtheta, par.phaseHamiltonian, par.lfAlpha);
dduNew = (op.Ltheta * uNew(:)).';

% r_if 的排列顺序需要与 W(:) 一致：先所有 x，再 theta。
rBase = Kx(:) - rho(:);
rIfMat = bsxfun(@plus, rBase, -HhatOld(:).' - eps*dduNew(:).');
rIf = rIfMat(:);

A0 = eps*speye(N) - dt*DblockLx - dt*eps^2*Lth + dt*2*eps*T;

switch lower(reactionMode)
    case {'patankar','production-destruction','pd','positive','positivity'}
        % production-destruction 型处理：正生产用 qW，负消耗用 W^{n+1}。
        % 这比把正 r_if 全隐式放进对角线更稳健。
        rPlus = max(rIf, 0);
        rMinus = max(-rIf, 0);
        A = A0 + dt*spdiags(rMinus,0,N,N);
        b = eps*qW(:) + dt*(rPlus .* qW(:));
    case {'implicit','linear-implicit','fully-implicit'}
        % 旧的线性隐式形式。小 eps 且 dt 较大时可能失去 M-matrix 结构。
        A = A0 - dt*spdiags(rIf,0,N,N);
        b = eps*qW(:);
    otherwise
        error('Unknown reactionDiscretization "%s".', reactionMode);
end

%--------------------------------------------------------------
% 3. 可选矩阵病态监测。condest 对大矩阵会增加耗时，因此默认关闭。
%--------------------------------------------------------------
info.matrixCondEst = NaN;
info.matrixRcondEst = NaN;
if local_get_logical(par, 'monitorMatrixCondition', false)
    try
        cA = condest(A);
        info.matrixCondEst = cA;
        info.matrixRcondEst = 1/cA;
    catch
        info.matrixCondEst = NaN;
        info.matrixRcondEst = NaN;
    end
end

Wvec = A \ b;
Wnew = reshape(real(Wvec), nx, Ktheta);

% 极小负值通常来自线性代数舍入；这里截断并记录。
tolNeg = local_get_num(par, 'negativeClipTol', 1e-13);
negMask = Wnew < 0;
info.numNegativeW = nnz(negMask);
info.minWBeforeClip = min(Wnew(:));
if any(negMask(:))
    smallNeg = negMask & (Wnew > -tolNeg*max(1,max(abs(Wnew(:)))));
    Wnew(smallNeg) = 0;
end
info.rIfMax = max(rIf);
info.rIfMin = min(rIf);
info.rIfPlusMax = max(max(rIf),0);
info.rIfAbsMax = max(abs(rIf));
info.WmaxAfter = max(Wnew(:));
info.WminAfter = min(Wnew(:));
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end

function val = local_get_num(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name)) && isscalar(s.(name))
    val = double(s.(name));
end
end

function tf = local_get_logical(s, name, defaultValue)
tf = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    tf = logical(s.(name));
end
end
