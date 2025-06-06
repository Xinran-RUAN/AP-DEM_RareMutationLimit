function u_new = solve_u(u, u_temp, B, H, d_theta, dt)
u_plus = [u_temp(2:end), u_temp(1)]; % ����theta�����ڱ߽�
u_minus = [u_temp(end), u_temp(1:end-1)];
% �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
pL = (u_temp - u_minus) / d_theta; % ���� u_j - u_{j-1}
pR = (u_plus - u_temp) / d_theta; % �Ҳ�� u_{j+1} - u_j
grad_sq = min(pL.^2, pR.^2) .* (pL .* pR >= 0) + ...
          pL .* 0 .* (pL .* pR < 0);
grad_sq(1) = pR(1).^2;
%    grad_sq1 = (Du(u) / d_theta).^2;
%     pp = spline([theta, 1], [u, u(1)]);
%     pp_der = fnder(pp, 1);
%     du = ppval(pp_der, theta);
%     grad_sq = du.^2;
% ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler
u_new = B\(u + dt * (grad_sq - H))';

u_new = u_new';
end

