function nNew = update_direct_density_fd_1d(n, D, Kx, eps, dt, op, par)
%UPDATE_DIRECT_DENSITY_FD_1D One step for the original density equation.
%   The default reaction treatment is controlled by par.reactionDiscretization:
%     'patankar' : positivity-preserving production-destruction reaction;
%     'implicit' : old linear implicit reaction n^{n+1}(K-rho^n).
%   If par is not supplied, the old implicit discretization is used.
Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;

if nargin < 7 || isempty(par)
    par = struct('reactionDiscretization','implicit');
end
if isfield(par, 'reactionDiscretization') && ~isempty(par.reactionDiscretization)
    reactionMode = char(par.reactionDiscretization);
else
    reactionMode = 'implicit';
end

rho = op.dtheta * sum(n, 2);
DblockLx = kron(spdiags(D(:),0,Ktheta,Ktheta), op.Lx);
Lth = op.Ltheta_big;

rBase = Kx(:) - rho(:);
rBaseVec = repmat(rBase, Ktheta, 1);

switch lower(reactionMode)
    case {'logistic-exact','exact-logistic','exact','if-exact'}
        % 可选 direct solver 的兼容处理：先对局部 logistic 反应作精确乘子，
        % 再隐式处理扩散。WKB 主求解器不使用这里的 direct n 更新。
        [logG, ~] = src.reaction_log_multiplier_1d(Kx, rho, eps, dt, par);
        nIf = bsxfun(@times, n, exp(min(max(logG(:),-745),700)));
        A = eps*speye(N) - dt*DblockLx - dt*eps^2*Lth;
        b = eps*nIf(:);
    case {'patankar','production-destruction','pd','positive','positivity'}
        rPlusVec = max(rBaseVec, 0);
        rMinusVec = max(-rBaseVec, 0);
        A = eps*speye(N) - dt*DblockLx - dt*eps^2*Lth + dt*spdiags(rMinusVec,0,N,N);
        b = eps*n(:) + dt*(rPlusVec .* n(:));
    case {'implicit','linear-implicit','fully-implicit'}
        A = eps*speye(N) - dt*DblockLx - dt*eps^2*Lth - dt*spdiags(rBaseVec,0,N,N);
        b = eps*n(:);
    otherwise
        error('Unknown reactionDiscretization "%s".', reactionMode);
end

nNew = reshape(A \ b, nx, Ktheta);
nNew = real(nNew);
end
