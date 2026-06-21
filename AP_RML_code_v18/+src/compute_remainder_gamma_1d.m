function gammaR = compute_remainder_gamma_1d(W, u, D, eps, op, par, variant)
%COMPUTE_REMAINDER_GAMMA_1D 加权正 remainder 诊断量。
%   使用 log-stabilized 重构，因此不要求 max(u)=0。
if nargin < 7 || isempty(variant)
    variant = 'split';
end
[~, n] = src.reconstruct_rho_1d(W, u, eps, op, 'direct-log');
Ktheta = op.Ntheta;
nx = op.nx;

switch lower(variant)
    case {'density-compatible','dc','density-compatible-if','dc-if','density-if','integrating-factor','density-compatible-semidiscrete','dc-semidiscrete','density-compatible-old','dc-old'}
        Rvec = eps^2 * (op.Ltheta_big * n(:));
        R = reshape(Rvec, nx, Ktheta);
    otherwise
        Hhat = src.numerical_hamiltonian(u, op.dtheta, par.phaseHamiltonian, par.lfAlpha);
        ddu = (op.Ltheta * u(:)).';
        ddW = reshape(op.Ltheta_big * W(:), nx, Ktheta);
        T = src.upwind_transport_matrix(u, op);
        DthetaF = reshape(T * W(:), nx, Ktheta);
        % E = exp(u/eps) 可能因 gauge 选择而上溢；用 E = n/W 的稳定形式。
        Eeff = n ./ max(W, realmin);
        R = eps^2*ddW.*Eeff - eps*bsxfun(@times, n, ddu) ...
            - bsxfun(@times, n, Hhat) - 2*eps*DthetaF.*Eeff;
end
RD = op.dtheta * sum(bsxfun(@rdivide, R, D), 2);
m = op.dtheta * sum(n, 2);
gammaR = max(max(real(RD),0) ./ (m + 1e-14));
end
