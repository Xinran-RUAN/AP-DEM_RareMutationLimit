function [Wnew, info] = update_w_density_if_fd_1d(W, uOld, uNew, D, Kx, rho, eps, dt, op, par)
%UPDATE_W_DENSITY_IF_FD_1D Time-density-compatible WKB amplitude update.
%   LEGACY comparison routine. This matrix form is not the small-epsilon default motivated by the
%   identity n = W*exp(u/eps).  Let
%
%       q_k = exp((uOld_k-uNew_k)/eps) = E^n_k/E^{n+1}_k.
%
%   The default update discretizes, after multiplication by E^{n+1}, the
%   density-level equation
%
%       eps (n^{n+1}-n^n)/dt - D L_x n^{n+1} - eps^2 L_theta n^{n+1}
%          = (K-rho^n) n^{n+1}
%
%   or, with par.reactionDiscretization='patankar', the positivity-preserving
%   production-destruction form
%
%       rhs = (K-rho^n)_+ n^n - (K-rho^n)_- n^{n+1}.
%
%   This removes the time-direction WKB chain-rule defect caused by separately
%   time-discretizing u and W.
info = struct();
Ktheta = op.Ntheta;
nx = op.nx;
N = nx*Ktheta;

if nargin < 10 || isempty(par)
    par = struct();
end
reactionMode = local_get_string(par, 'reactionDiscretization', 'patankar');
expClip = local_get_num(par, 'expClip', 700);

% x diffusion: each theta block uses D(theta_k).
DblockLx = kron(spdiags(D(:),0,Ktheta,Ktheta), op.Lx);

% theta diffusion acting on the reconstructed density, divided by E^{n+1}.
Btheta = src.density_compatible_theta_matrix(uNew, eps, op);

% qW is the old density divided by E^{n+1}; i.e. q_k W^n_{j,k}.
q = src.safe_exp((uOld(:).' - uNew(:).')/eps, expClip);
qW = bsxfun(@times, W, q);

rBase = Kx(:) - rho(:);
rBaseVec = repmat(rBase, Ktheta, 1);

switch lower(reactionMode)
    case {'patankar','production-destruction','pd','positive','positivity'}
        rPlusVec = max(rBaseVec, 0);
        rMinusVec = max(-rBaseVec, 0);
        A = eps*speye(N) - dt*DblockLx - dt*eps^2*Btheta + dt*spdiags(rMinusVec,0,N,N);
        b = eps*qW(:) + dt*(rPlusVec .* qW(:));
    case {'implicit','linear-implicit','fully-implicit'}
        A = eps*speye(N) - dt*DblockLx - dt*eps^2*Btheta - dt*spdiags(rBaseVec,0,N,N);
        b = eps*qW(:);
    otherwise
        error('Unknown reactionDiscretization "%s".', reactionMode);
end

Wnew = reshape(A \ b, nx, Ktheta);
Wnew = real(Wnew);
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
