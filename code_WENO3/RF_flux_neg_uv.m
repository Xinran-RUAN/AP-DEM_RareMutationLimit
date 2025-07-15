function Hhat = RF_flux_neg_uv(u_plus, u_minus, v_plus, v_minus)
% Roe RF flux for H(u,v) = -u v (�ϸ񰴹�ʽ)
H = @(u,v) -u .* v;

% ƽ��ֵ���� upwind �ж�
u_avg = 0.5 * (u_plus + u_minus);
v_avg = 0.5 * (v_plus + v_minus);

% �ж� sign change
H1_nochange = v_plus .* v_minus >= 0;
H2_nochange = u_plus .* u_minus >= 0;

% upwind u*��v*
u_star = u_plus;
u_star(v_avg < 0) = u_minus(v_avg < 0);

v_star = v_plus;
v_star(u_avg < 0) = v_minus(u_avg < 0);

% ��ʼ��
Hhat = zeros(size(u_plus));

% �ֶι�ʽ
mask1 = H1_nochange & H2_nochange;
Hhat(mask1) = H(u_star(mask1), v_star(mask1));

mask2 = ~H1_nochange & H2_nochange;
alpha = 0; % ��Ϊ H11 = 0
Hhat(mask2) = H(0.5*(u_plus(mask2)+u_minus(mask2)), v_star(mask2)) ...
              - alpha * (u_plus(mask2)-u_minus(mask2)) .* (u_plus(mask2)-u_minus(mask2)) / 2;

mask3 = H1_nochange & ~H2_nochange;
beta = 0; % ��Ϊ H22 = 0
Hhat(mask3) = H(u_star(mask3), 0.5*(v_plus(mask3)+v_minus(mask3))) ...
              - beta * (v_plus(mask3)-v_minus(mask3)) .* (v_plus(mask3)-v_minus(mask3)) / 2;

mask4 = ~(mask1 | mask2 | mask3);
% LLF ��д��ƽ��
Hhat(mask4) = 0.5 * (H(u_plus(mask4), v_plus(mask4)) + H(u_minus(mask4), v_minus(mask4)));
end
