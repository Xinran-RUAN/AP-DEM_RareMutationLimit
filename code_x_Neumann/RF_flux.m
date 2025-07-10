function Hhat = RF_flux(u_plus, u_minus)
% ʸ���� Roe entropy fix flux for H(u) = -u^2

% ���� H �� H'
H = @(u) -u.^2;
Hp = @(u) -2*u;

% Roe ƽ��
ubar = 0.5 * (u_plus + u_minus);

% ����˵�
pmin = min(u_minus, u_plus);
pmax = max(u_minus, u_plus);

Hpmin = Hp(pmin);
Hpmax = Hp(pmax);

% upwind ѡ��
u_star = u_plus;
u_star(Hp(ubar) >= 0) = u_minus(Hp(ubar) >= 0);

% �ж� sign change
no_sign_change = (Hpmin >= 0 & Hpmax >= 0) | (Hpmin <= 0 & Hpmax <= 0);

% ��ʼ��
Hhat = zeros(size(u_plus));

% �� sign change
Hhat(no_sign_change) = H(u_star(no_sign_change));

% �� sign change
alpha = 1;
du = u_plus - u_minus;
Hhat(~no_sign_change) = H(0.5*(u_plus(~no_sign_change)+u_minus(~no_sign_change))) ...
    - alpha * du(~no_sign_change) .* du(~no_sign_change) / 2;
end



