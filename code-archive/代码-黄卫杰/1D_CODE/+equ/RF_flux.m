function Hhat = RF_flux(u_plus, u_minus)
% 矢量化 Roe entropy fix flux for H(u) = -u^2

% 定义 H 和 H'
H = @(u) -u.^2;
Hp = @(u) -2*u;

% Roe 平均
ubar = 0.5 * (u_plus + u_minus);

% 区间端点
pmin = min(u_minus, u_plus);
pmax = max(u_minus, u_plus);

Hpmin = Hp(pmin);
Hpmax = Hp(pmax);

% upwind 选择
u_star = u_plus;
u_star(Hp(ubar) >= 0) = u_minus(Hp(ubar) >= 0);

% 判断 sign change
no_sign_change = (Hpmin >= 0 & Hpmax >= 0) | (Hpmin <= 0 & Hpmax <= 0);

% 初始化
Hhat = zeros(size(u_plus));

% 无 sign change
Hhat(no_sign_change) = H(u_star(no_sign_change));

% 有 sign change
alpha = 2 * max(abs(u_plus), abs(u_minus));
du = u_plus - u_minus;
Hhat(~no_sign_change) = H(0.5*(u_plus(~no_sign_change)+u_minus(~no_sign_change))) ...
    - alpha(~no_sign_change) .* du(~no_sign_change) / 2;
end



