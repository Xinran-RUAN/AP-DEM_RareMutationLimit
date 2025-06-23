function u_new = solve_u(u, u_temp, B, H, d_theta, dt, theta, ol)
u_plus = [u_temp(2:end), u_temp(1)]; % ����theta�����ڱ߽�
u_minus = [u_temp(end), u_temp(1:end-1)];
% �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
if ol == 0
    %% ȡ����ֵС�ĵ���
    pL = (u_temp - u_minus) / d_theta; % ���� u_j - u_{j-1}
    pR = (u_plus - u_temp) / d_theta; % �Ҳ�� u_{j+1} - u_j
    grad_sq = min(pL.^2, pR.^2) .* (pL .* pR >= 0) + ...
        pL .* 0 .* (pL .* pR < 0);
else
    %% spline��ֵ����
    du = 0 .* u;
    pp = spline(theta(2:end), u(2:end));
    pp_der = fnder(pp, 1);
    du(2:end) = ppval(pp_der, theta(2:end));
    grad_sq = du.^2;
    pR = (u_plus - u_temp) / d_theta; % �Ҳ�� u_{j+1} - u_j
    du(1) = pR(1);
end
% ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler
u_new = B\(u + dt * (grad_sq - H))';

u_new = u_new';
end

