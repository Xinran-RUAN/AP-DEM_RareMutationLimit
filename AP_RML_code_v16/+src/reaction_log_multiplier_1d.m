function [logG, info] = reaction_log_multiplier_1d(Kx, rho, eps, dt, par)
%REACTION_LOG_MULTIPLIER_1D WKB 振幅中的反应积分因子。
%   该函数返回每个 x 网格点上的 logG，使得零维反应快层近似写成
%       W <- G(x) W,   G=exp(logG).
%   这是对 W 的乘法更新，仍然处在 WKB 框架内；并没有把主求解器改成 n 方程。
%
%   mode='logistic-exact'：对局部总密度 ODE
%       eps rho_t = (K-rho) rho
%   使用精确解给出 rho^{new}/rho^old，再把同一乘子作用到所有 theta。
%   这样可去掉正增长项导致的 dt/eps 保正限制。

mode = local_get_string(par, 'reactionDiscretization', 'logistic-exact');
K = real(Kx(:));
rho0 = max(real(rho(:)), realmin);

tinyK = 1e-14;
info = struct();
info.mode = mode;

switch lower(mode)
    case {'logistic-exact','exact-logistic','exact','if-exact'}
        G = ones(size(rho0));
        idx = K > tinyK;
        if any(idx)
            a = (K(idx)./rho0(idx) - 1) .* exp(-K(idx)*dt/eps);
            den = 1 + a;
            den = max(den, realmin);
            rho1 = K(idx) ./ den;
            G(idx) = rho1 ./ rho0(idx);
        end
        if any(~idx)
            % K=0 时为 eps rho_t = -rho^2。
            G(~idx) = 1 ./ (1 + dt*rho0(~idx)/eps);
        end
        logG = log(max(G, realmin));
    case {'patankar','production-destruction','pd','positive','positivity'}
        r = K - rho0;
        rp = max(r,0);
        rm = max(-r,0);
        logG = log1p(dt*rp/eps) - log1p(dt*rm/eps);
    case {'implicit','linear-implicit','fully-implicit'}
        r = K - rho0;
        den = 1 - dt*r/eps;
        % 线性隐式正增长在 den<=0 时不保正，这里只做保护并记录。
        bad = den <= realmin;
        den(bad) = realmin;
        logG = -log(den);
        info.numBadImplicitDenominator = nnz(bad);
    case {'none','off'}
        logG = zeros(size(rho0));
    otherwise
        error('Unknown reactionDiscretization "%s".', mode);
end

logG(~isfinite(logG)) = 0;
info.logGMax = max(logG);
info.logGMin = min(logG);
info.GMax = exp(min(info.logGMax, 700));
info.GMin = exp(max(info.logGMin, -745));
end

function val = local_get_string(s, name, defaultValue)
val = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = char(s.(name));
end
end
