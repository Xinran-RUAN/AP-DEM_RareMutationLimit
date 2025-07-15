function u_new = solve_u(u, u_temp, B, H, d_theta, dt, theta, ol)
%% �����ʽȫ������������u��ʾ��һ��ʱ�䲽��ֵ��u_temp��ʾ���������е�ֵ
u_plus = [u_temp(2:end), u_temp(1)]; 
u_minus = [u_temp(end), u_temp(1:end-1)];
% �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��

%% spline��ֵ����
pp = spline(theta(2:end-1), u(2:end-1));
pp_der = fnder(pp, 1);
du = ppval(pp_der, theta(2:end-1));
grad_sq = du.^2;
   
% ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler
u_new = B\(u(2:end-1) + dt * (grad_sq - H(2:end)))';

u_new = [u_new(1); u_new; u_new(end)]';
end

